# -*- coding: UTF-8 -*-
import os
import glob
import re
import threading
import shutil

class HMM(object):

    def __init__(self, nipt):
        self.indir = nipt.pca_dir
        self.outdir = nipt.hmm_dir
        self.zscore = nipt.hmm_zscore
        self.length = nipt.hmm_length
        self.rate = nipt.hmm_rate
        self.zsd = nipt.hmm_zsd
        self.xhmm_exe = nipt.hmm_exe
        self.disease_file = nipt.hmm_disease
        self.flowcell = nipt.flowcell
        self.fid = self.flowcell.split('_')[-1]
        self.chrs = nipt.chrs
        self.samples = nipt.sample_deal
        self.param_file = "/".join([self.indir, 'param.txt'])
        self.report_file_prefix = "/".join([self.outdir,"CNV_report"])
        self._disease = None
        self._disease_names = None
        self._rd_files = None
        self._fd_matrix = None
        self._xcnv_files = None
        self._cnv_report_detail_head = None
        self._cnv_result = None
        self._hmm_threads = None
        self._simplify_cnv_result = None
        self._cnv_report_txt_head = None
        self._disease_report = None

    def get_disease_report(self):
        if self._disease_report:
            return self._disease_report
        disease_sample = {}
        for sample, cnvreports in self.simplify_cnv_result.iteritems():
            disease_sample[sample] = {}
            for chrnum, cnv in cnvreports.iteritems():
                branch = cnv.split("\n")
                branch.pop() # 去除""
                for bra in branch:
                    items = re.compile("\s+").split(bra)
                    dis = items[-1].split(";")
                    for muti in dis:
                        if items[-2] == 'mat':
                            disease_sample[sample] = { muti : {'mat' : 'positive' } }
                        else:
                            disease_sample[sample] = { muti : {'fetal' : 'positive' } }
        self._disease_report = disease_sample
        return disease_sample
    disease_report = property(get_disease_report)

    def get_simplify_cnv_result(self):
        if self._simplify_cnv_result:
            return self._simplify_cnv_result
        inputfile = ".".join([self.report_file_prefix, 'detail'])
        result = {}
        with open(inputfile, 'r') as f:
            for line in f:
                if re.compile("^\s*$").search(line):
                    continue
                items = re.compile("\s+").split(line.strip())
                if line.startswith("SAMPLE"):
                    if not self._cnv_report_txt_head:
                        cnv_report_txt_head = "\t".join([items[0], items[1], items[2], items[3], items[4],
                            items[13], items[15], items[16], items[17]])
                        self._cnv_report_txt_head = "".join([cnv_report_txt_head, "\n"])
                else:
                    if float(items[3]) < 2000:
                        continue
                    else:
                        cnv_content = "\t".join([items[0], items[1], items[2], items[3], items[4],
                            items[13], items[15], items[16], items[17]])
                        simplify_cnv_result_sample_chr_single = "".join([cnv_content, "\n"])
                        if result.has_key(items[0]):
                            if result[items[0]].has_key(items[4]):
                                result[items[0]][items[4]] = "".join([result[items[0]][items[4]], simplify_cnv_result_sample_chr_single])
                            else:
                                result[items[0]][items[4]] = simplify_cnv_result_sample_chr_single
                        else:
                            result[items[0]] = {}
                            result[items[0]][items[4]] = simplify_cnv_result_sample_chr_single
            self._simplify_cnv_result = result
            return result
    simplify_cnv_result = property(get_simplify_cnv_result)

    def get_cnv_report_txt_head(self):
        if self._cnv_report_txt_head:
            return self._cnv_report_txt_head
        self.get_simplify_cnv_result()
        return self._cnv_report_txt_head
    cnv_report_txt_head = property(get_cnv_report_txt_head)

    def get_hmm_threads(self):
        if self._hmm_threads:
            return self._hmm_threads
        threads = []
        for file_get in self.rd_files:
            filename = os.path.basename(file_get)
            prefix = filename.replace('.normalized.zscore.matrix','')
            raw_rd = file_get.replace('normalized.zscore.matrix','redis_scale.txt')
            if not os.path.isfile(raw_rd):
                raw_rd = raw_rd.replace('redis_scale.txt', 'filter_na.matrix')
                command = [
                    self.xhmm_exe, '--discover -p', self.param_file,
                    '-r', file_get, '-R', raw_rd,
                    '-c', '/'.join([self.outdir, '.'.join([prefix,'xcnv'])]),
                    '-a', '/'.join([self.outdir, '.'.join([prefix,'aux_xcnv'])]),
                    '-s', '/'.join([self.outdir, prefix]),
                    '>/dev/null 2>&1'
                ]
                print ' '.join(command)
                threads.append(threading.Thread(target=os.system,args=(' '.join(command),)))
        self._hmm_threads = threads
        return threads
    hmm_threads = property(get_hmm_threads)

    def get_cnv_result(self):
        if self._cnv_result:
            return self._cnv_result
        cnv_result = {}
        for xcnv in self.xcnv_files:
            with open(xcnv, 'r') as f:
                for line in f:
                    if line.startswith("SAMPLE"):
                        if not self._cnv_report_detail_head:
                            self._cnv_report_detail_head = "\t".join([line.strip(),'cffDNA%','ORIGIN',"DISEASE\n"])
                    else:
                        items = line.split("\t")
                        items[1] = 'loss' if items[1] == 'DEL' else 'gain'
                        if float(items[3]) < 1000:
                            continue
                        if cnv_result.has_key(items[0]):
                            cnv_result[items[0]][items[4]] = "".join([cnv_result[items[0]][items[4]], "\t".join(items)])\
                                    if cnv_result[items[0]].has_key(items[4]) else "\t".join(items)
                        else:
                            cnv_result[items[0]] = {}
                            cnv_result[items[0]][items[4]] = "\t".join(items)
        self._cnv_result = cnv_result
        return cnv_result
    cnv_result = property(get_cnv_result)

    def get_cnv_report_detail_head(self):
        if self._cnv_report_detail_head:
            return self._cnv_report_detail_head
        self.get_cnv_result()
        return self._cnv_report_detail_head
    cnv_report_detail_head = property(get_cnv_report_detail_head)

    def get_xcnv_files(self):
        if self._xcnv_files:
            return self._xcnv_files
        xcnv_files = glob.glob("/".join([self.outdir,'*.xcnv']))
        xcnv_files.sort()
        self._xcnv_files = xcnv_files
        return xcnv_files
    xcnv_files = property(get_xcnv_files)

    def get_fd_matrix(self):
        if self._fd_matrix:
            return self._fd_matrix
        result = {}
        for filename in self.rd_files:
            fd_matrix_file = filename.replace("normalized.zscore.matrix","fd.matrix")
            chr_num = re.compile("chr(\d+)").search(os.path.basename(filename)).group(1)
            result[chr_num] = fd_matrix_file
        self._fd_matrix = result
        return result
    fd_matrix = property(get_fd_matrix)

    def parse_disease(self):
        if self._disease:
            return self._disease
        result = {}
        disease_name = []
        with open(self.disease_file, 'r') as f:
            for line in f:
                line = line.strip()
                if line.startswith("#") or re.match("^\s+$", line):
                    continue
                items = re.compile("\s+").split(line)
                name = items.pop(0)
                disease_name.append(name)
                result[name] = items
        result_sorted = sorted(result.iteritems(), key=lambda d:d[0]);
        self._disease = result_sorted
        if not self._disease_names:
            self._disease_names = disease_name
        return result_sorted
    disease = property(parse_disease)

    def get_disease_names(self):
        if self._disease_names:
            return self._disease_names
        self.parse_disease()
        return self._disease_names
    disease_names = property(get_disease_names)

    def get_rd_files(self):
        if self._rd_files:
            return self._rd_files
        result = []
        for chr_num in self.chrs:
            file_get = glob.glob("/".join([self.indir, ''.join(["*.chr",str(chr_num),'.RD.normalized.zscore.matrix'])]))[0]
            result.append(file_get)
        self._rd_files = result
        return result
    rd_files = property(get_rd_files)

    def __str__(self):
        return str(vars(self))

    def generate_param_file(self):
        with open(self.param_file, 'w') as f:
            content = "\t".join([
                str(self.rate), str(self.length), '70',
                str(-self.zscore), str(self.zsd), '0',
                str(self.zsd), str(self.zscore), str(self.zsd)
            ])
            f.write(content)

    def hmm_analyse(self):
        print "HMM分析开始"
        self.generate_param_file()
        # check outdir
        if (os.path.exists(self.outdir)):
            pass
            old_staff = glob.glob('/'.join([self.outdir,'*']))
            for staff in old_staff:
                os.remove(staff) if os.path.isfile(staff) else shutil.rmtree(staff)
        else:
            os.makedirs(self.outdir)
        for t in self.hmm_threads:
            t.setDaemon(True)
            t.start()
        for t in self.hmm_threads:
            t.join()
        print "调用XHMM步骤完毕"
        # merge the xhmm result by chromosome, this is a detail result
        # the cnv region length more than 1M
        # 生成CNV_report.detail
        self.generate_cnv_report_detail()
        # 生成CNV_report.txt
        self.generate_cnv_report_txt()
        # 生成CNV_report.ForAnno -- 暂缓
        # 生成CNV_report.seven.dis.csv
        self.generate_cnv_report_seven_dis_csv()

    def generate_cnv_report_seven_dis_csv(self):
        other_mat_disorder = 'Fetal chromosomal aneuploidy (>= 2M bases)'
        other_fetal_disorder = 'Maternal chromosomal aneuploidy (>= 2M bases)'
        cnv_report_seven_dis_csv = ".".join([self.report_file_prefix,'seven.dis.csv'])
        disease_sample = self.disease_report
        with open(cnv_report_seven_dis_csv,'w') as f:
            header = ["".join([i,"(fetal/maternal)"]) for i in self.disease_names]
            header.insert(0,"SAMPLE")
            header.append(other_mat_disorder)
            header.append(other_fetal_disorder)
            f.write(','.join(header))
            f.write("\n")
            for sample in sorted(self.samples):
                print sample, " -------> generate Cnv_report.seven.dis.csv"
                sample_dis_info = [sample]
                for dis in self.disease_names:
                    if disease_sample.has_key(sample) and disease_sample[sample].has_key(dis):
                        mat = 'positive' if disease_sample[sample][dis].has_key('mat') else 'negative'
                        fetal = 'posivite' if disease_sample[sample][dis].has_key('fetal') else 'negative'
                        sample_dis_info.append("/".join([mat, fetal]))
                    else:
                        sample_dis_info.append("/".join(['negative', 'negative']))
                if disease_sample.has_key(sample) and disease_sample[sample].has_key('other'):
                    mat_other = "positive" if disease_sample[sample]['other'].has_key('mat') else 'negative'
                    fetal_other = "positive" if disease_sample[sample]['other'].has_key('fetal') else 'negative'
                    mat_other = mat_other + "\n"
                    sample_dis_info.append(fetal_other)
                    sample_dis_info.append(mat_other)
                else:
                    sample_dis_info.append("negative")
                    sample_dis_info.append("negative\n")
                f.write(",".join(sample_dis_info))

    def generate_cnv_report_txt(self):
        cnv_report_txt = ".".join([self.report_file_prefix, "txt"])
        with open(cnv_report_txt,'w') as f:
            f.write(self.cnv_report_txt_head)
            for sample, cnvreports in sorted(self.simplify_cnv_result.iteritems(), key=lambda d: d[0]):
                print sample, ' -------> generate Cnv_report.txt'
                for chrnum, cnvreport in sorted(cnvreports.iteritems(), key=lambda k: k[0]):
                    f.write(cnvreport)

    def generate_cnv_report_detail(self):
        cnv_report_detail = ".".join([self.report_file_prefix, 'detail'])
        #cnv_result_sorted = sorted(self.cnv_result.iteritems(), key=lambda d:d[0]);
        #print cnv_result_sorted
        with open(cnv_report_detail, 'w') as f:
            f.write(self.cnv_report_detail_head)
            for (sample, cnvs) in sorted(self.cnv_result.iteritems(), key=lambda d:d[0]):
                print sample, " -------> generate Cnv_report.detail"
                for (chr_num, cnv) in sorted(cnvs.iteritems(), key=lambda d:d[0]):
                    #print chr_num, '#--->', cnv
                    #print self.fd_matrix[chr_num]
                    # count_cffdna
                    cnv_add_cffdna = self.count_cffdna(sample, chr_num, cnv)
                    #print cnv_add_cffdna
                    total_record_single = self.get_seven_dis_info(cnv_add_cffdna);
                    #print total_record_single
                    f.write(total_record_single)


    def get_seven_dis_info(self, cnv_add_cffdna):
        cnvs = cnv_add_cffdna.split("\n")
        result = []
        for cnv in cnvs:
            if cnv == '':
                continue
            items = re.compile("\s+").split(cnv)
            (chrnum, start, end) = map(int,re.compile(":|-").split(items[2]))
            dis_count = []
            for dis_name, dis_content in self.disease:
                #print dis_name, "#------->", dis_content
                if chrnum == int(dis_content[0])\
                and items[1] == dis_content[3]\
                and start < int(dis_content[2])\
                and end > int(dis_content[1]):
                    n_start = start if start > int(dis_content[1]) else int(dis_content[1])
                    n_end = end if end < int(dis_content[2]) else int(dis_content[2])
                    if (end-start) >= 2000000:
                        dis_count.append(dis_name)
            if len(dis_count) < 1:
                result.append("\t".join([cnv, "other"]))
            else:
                result.append("\t".join([cnv, ";".join(dis_count)]))
        result.append("")
        return "\n".join(result)

    def count_cffdna(self, sample, chrnum, cnv):
        rawfile = self.fd_matrix[chrnum]
        use_regions = None
        use_rc = None
        with open(rawfile, 'r') as f:
            use_regions = re.compile("\s+").split(re.sub("^\s+","", f.readline().strip()))
            for line in f:
                if line.startswith(sample):
                    items = re.compile("\s+").split(line.strip())
                    items.pop(0)
                    use_rc = items
                    break
        records = cnv.split("\n")
        records.pop() # remove last item: ""
        #print records
        cnv_nu = {}
        cnv_sum = {}
        normal_nu = 0
        normal_sum = 0
        for region, rc in zip(use_regions, use_rc):
            #print region, "#----->", rc
            stat = 0
            (chrnum, start, end) = map(int,re.compile(":|-").split(region))
            rc = abs(float(rc))
            for record in records:
                items = re.compile("\s+").split(record)
                (n_chr, n_start, n_end) = map(int, re.compile(":|-").split(items[2]))
                if n_start <= end and n_end >= start:
                    stat += 1
                    cnv_nu[items[2]] = cnv_nu[items[2]] + 1 if cnv_nu.has_key(items[2]) else 1
                    cnv_sum[items[2]] = cnv_sum[items[2]] + rc if cnv_sum.has_key(items[2]) else rc
            if stat == 0:
                normal_nu += 1
                normal_sum += rc
        a_cffdna = []
        for record  in records:
            items = re.compile("\s+").split(record)
            cffdna = abs(2- (cnv_sum[items[2]] / cnv_nu[items[2]]) / (normal_sum / normal_nu) * 2)
            origin = 'fetal' if abs(float(items[13])) < 2.5 and cffdna < 0.6 else 'mat'
            a_cffdna.append("\t".join([record, "%.2f" % (cffdna*100), origin]))
        a_cffdna.append("")
        return "\n".join(a_cffdna)


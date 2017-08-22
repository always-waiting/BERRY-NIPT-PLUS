# -*- coding: UTF-8 -*-
import os
import glob
import re
import shutil
from multiprocessing.dummy import Pool as ThreadPool

class PCA(object):
    """PCA分析类"""
    def __init__(self, niptobj):
        self.indir = niptobj.trim_dir
        self.outdir = niptobj.pca_dir
        self.flowcell = niptobj.flowcell
        self.fid = self.flowcell.split('_')[-1]
        self.chrs = niptobj.chrs
        self.samples = niptobj.sample_deal
        self.rm_mean_pve = niptobj.rm_mean_pve if niptobj.rm_mean_pve else 0.95
        self.rmpcs = niptobj.rmPCs
        self.redisseq = niptobj.redisseq
        self.param_dir = niptobj.pca_param_dir
        self.merge = niptobj.merge_block
        self.parallel = niptobj.parallel

    def make_rd_matrix(self):
        print "RD矩阵生成"

        if not os.path.exists(self.indir):
            raise Exception(''.join([self.indir, ' not Exist!']))

        if not os.path.exists(self.outdir):
            os.makedirs(self.outdir)
        else:
            old_staff = glob.glob("/".join([self.outdir,"*"]))
            for staff in old_staff:
                os.remove(staff) if os.path.isfile(staff) else shutil.rmtree(staff)

        #print vars(self)
        samples = list(self.samples)
        samples.sort()
        all_samples = {}
        for num in self.chrs:
            all_samples[num] = []
        for sample in samples:
            filename = glob.glob('/'.join([self.indir,''.join([sample,'*.filter.*.txt'])]))
            if len(filename) == 0:
                raise Exception(''.join([sample, '.filter.*.txt is not found']))
            #print filename
            with open(filename[0],'r') as f:
                header = f.readline()
                total = int(re.search('\d+', header).group(0))
                #print total
                for line in f:
                    lines = line.split('\t')
                    chr_num = int(re.search('\d+', lines.pop(0)).group(0))
                    if chr_num in self.chrs:
                        #print ''.join([str(chr_num), ' should be counted'])
                        lines = map((lambda x: x if not x.isdigit() else "%.12f" % (float(x)*20000000/total)) ,lines)
                        lines.insert(0,sample)
                        all_samples[chr_num].append("\t".join(lines))
        #print all_samples[1]
        for key, value in all_samples.iteritems():
            content = "".join(value)
            outfile = "/".join([self.outdir, ''.join([self.fid,'.chr',str(key),'.RD.matrix'])])
            with open(outfile,'w') as f:
                f.write(content)
        print "RD矩阵生成完毕"

    def get_pca_r_option(self):
        r_option_list = ['--Merge', str(self.merge)]
        if self.redisseq:
            r_option_list.append('--redisseq')
            r_option_list.append(str(self.redisseq))
        if self.rmpcs:
            r_option_list.append('--rmPCs')
            r_option_list.append(str(self.rmpcs))
        else:
            r_option_list.append('--rm_mean_pve')
            r_option_list.append(str(self.rm_mean_pve))
        return ' '.join(r_option_list)

    pca_r_option = property(get_pca_r_option)

    def get_r_script(self):
        dirlocate = os.path.dirname(os.path.abspath(__file__))
        if self.param_dir:
            return '/'.join([dirlocate,'lib','R','pca_normalize_use_param.r'])
        else:
            return '/'.join([dirlocate,'lib','R','pca_normalize.r'])

    rscript = property(get_r_script)

    def pca(self):
        print "开始PCA分析"
        #print vars(self)
        print self.pca_r_option
        pool = ThreadPool(self.parallel)
        param_list = list(self.chrs)
        pool.map(self.pca_single, param_list)
        print "PCA分析完成"
        print "生成foldchange_sd.matrix"
        sdfiles = glob.glob('/'.join([self.outdir,'*.fd.mean_sd.matrix']))
        sdfiles.sort()
        foldchange_sd_content = []
        sample_dict = {}
        sample_list = []
        FD_SD_MAX_THRESHOLD = 0.08
        for filename in sdfiles:
            with open(filename, 'r') as f:
                if len(foldchange_sd_content) == 0:
                    items = f.readline().split('\t')
                    items.pop(0)
                    sample_list = items[:]
                    for sample in items:
                        sample_dict[sample] = []
                    items.insert(0,'Sample')
                    foldchange_sd_content.append('\t'.join(items))
                else:
                    f.readline()
                f.readline()
                items = f.readline().split('\t')
                items.pop(0)
                chr_str = re.search('chr\d+',f.name).group(0)
                chr_str = re.search('\d+', chr_str).group(0)
                fd_list = [ n for n,i in enumerate(items) if float(i)>FD_SD_MAX_THRESHOLD ]
                if len(fd_list) > 0:
                    for index in fd_list:
                        sample_dict[sample_list[index]].append(chr_str)
                items.insert(0,"".join([chr_str,"_sd"]))
                foldchange_sd_content.append('\t'.join(items))

        with open('/'.join([self.outdir, 'foldchange_sd_content.matrix']), 'w') as f:
            f.write(''.join(foldchange_sd_content))
        #print sample_dict
        FD_DISPERSE_NUMBER_THRESHOLD = 4
        sample_dict_sorted = sorted(sample_dict.iteritems(), key=lambda d:d[0]);
        with open('/'.join([self.outdir, 'constancy_sd.txt']), 'w') as f:
            content = ["Sample\tconstancy\tfreq\tcontent"]
            chr_total = len(self.chrs)
            for (key, value) in sample_dict_sorted:
                if len(value) >= FD_DISPERSE_NUMBER_THRESHOLD:
                    # 开始生成constancy_sd.txt
                    content.append("\t".join([key.strip(),'NO',"%d/%d" % (len(value),chr_total), ';'.join(value)]))
                else:
                    if len(value) == 0:
                        content.append("\t".join([key,'YES',"%d/%d" % (len(value), chr_total), '--']))
                    else:
                        content.append("\t".join([key,'YES',"%d/%d" % (len(value), chr_total), ';'.join(value)]))
            f.write("\n".join(content))



    def pca_single(self, chr_num):
        print chr_num
        rdfile = '/'.join([self.outdir,''.join([self.fid,'.chr',str(chr_num),'.RD.matrix'])])
        command_list = [\
            '/usr/bin/Rscript', self.rscript,\
            '--input', rdfile,\
            '--outdir', self.outdir,\
            '--Chr', str(chr_num),\
            self.pca_r_option\
        ]
        if self.param_dir:
            command_list.append('--prefile')
            refs = '/'.join([self.param_dir, ''.join(['hk.chr',str(chr_num),'.RD'])])
            command_list.append(refs)
        command = ' '.join(command_list)
        print command
        os.system(command)





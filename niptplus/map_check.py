# -*- coding: UTF-8 -*-
import os
import glob
import re
import shutil


class CheckData(object):
    def __init__(self, nipt):
        self.indir = nipt.indir
        self.outdir = nipt.map_dir
        self.samples = nipt.sample_deal
        self.min_uniq_reads = nipt.min_uniq_reads

    def check_map(self):
        if not os.path.exists(self.indir):
            print ' '.join([self.indir, 'not exists'])
            exit(1)
        if os.path.exists(self.outdir):
            old_staff = glob.glob(self.outdir)[0]
            shutil.rmtree(old_staff)
        lessdata = '/'.join([self.outdir, 'LessData'])
        os.makedirs(lessdata) # 上一步删除了oudir,所以必定没有LessData文件夹
        rd_files = glob.glob('/'.join([self.indir, '*.20K.txt']))
        gc_files = glob.glob('/'.join([self.indir, '*.20K.GC.txt']))
        pq_files = glob.glob('/'.join([self.indir, '*.PQ.txt']))
        def trim_files(data):
            trim_result =  filter(lambda x: re.split('\.|_',re.split('/', x)[-1])[0] in self.samples, data)
            trim_result.sort()
            return trim_result
        rd_files = trim_files(rd_files)
        gc_files = trim_files(gc_files)
        pq_files = trim_files(pq_files)
        if not len(rd_files) > 0\
        or (len(rd_files) != len(gc_files)\
        or len(gc_files) != len(pq_files)):
            print 'Sample files is wrong!'
            exit(1)
        with open('/'.join([self.outdir, 'mapping_result']), 'w') as mpf:
            mpf.write("\t".join(["Sample","Total_Rds","UniMap_Rds","Uniq%",
                "PCR_Dup","Redundancy%","UniMap_GC%","content","LessData\n"]))
            for filename in pq_files:
                pqfile = PQFile(filename, self.min_uniq_reads)
                pqfile.parse_file()
                pqfile.count_content()
                pqfile.count_lessdata()
                mpf.write(str(pqfile))
                if pqfile.lessdata == "YES":
                    shutil.copy(filename, lessdata)
                    shutil.copy(filename.replace("PQ.txt", "20K.GC.txt"), lessdata)
                    shutil.copy(filename.replace("PQ.txt", "20K.txt"), lessdata)
                else:
                    shutil.copy(filename, self.outdir)
                    shutil.copy(filename.replace("PQ.txt","20K.GC.txt"),self.outdir)
                    shutil.copy(filename.replace("PQ.txt","20K.txt"),self.outdir)


class PQFile(object):
    def __str__(self):
        return "\t".join([self.sample, str(self.total_rds), str(self.unimap_rds), self.uniq, str(self.pcr_dup), self.redundancy, self.unimap_gc, self.content, "".join([self.lessdata, "\n"])])

    def __init__(self, filename,  min_uniq_reads):
        self.filename = filename
        self.uniq = ''#
        self.unimap_gc = ''#
        self.unimap_rds = 0#
        self.total_rds = 0#
        self.pcr_dup = 0#
        self.redundancy = ''#
        self.sample = ''#
        self.content = ''
        self.lessdata = ''
        self.min_uniq_reads = min_uniq_reads

    def parse_file(self):
        print "解析", self.filename, "，存储结果"
        self.sample = self.filename.split("/")[-1].split(".")[0]
        with open(self.filename) as f:
            nonxym = 0
            total_gc = 0
            unimap_rds = 0
            totalmapped = 0
            for line in f:
                line = line.replace("NA","0")
                if line[0].isdigit():
                    items = line.strip().split()
                    if items[0] != 23 and items[0] != 24 and items[0] != 25:
                        nonxym += int(items[1])
                    total_gc += int(items[3])
                if line.startswith(">uniqMapped"):
                    unimap_rds = int(line.strip().split()[1])
                if line.startswith(">totalMapped"):
                    totalmapped = int(line.strip().split()[1])
                if line.startswith(">redundancy"):
                    self.redundancy = "%.3f" % (float(line.strip().split()[1])*100)
                if line.startswith(">no_N_Reads"):
                    self.total_rds = int(line.strip().split()[1])
            self.unimap_rds = unimap_rds
            self.pcr_dup = totalmapped - self.unimap_rds
            if self.total_rds != 0:
                self.uniq = "%.3f" % (100*self.unimap_rds/float(self.total_rds))
            else:
                self.uniq = 'NA'
            if self.unimap_rds != 0:
                self.unimap_gc = "%.3f" % (100*total_gc/float(36*self.unimap_rds))
            else:
                self.unimap_gc = "NA"

    def count_content(self):
        print "确定content"
        unimap_gc = float(self.unimap_gc)
        uniq = float(self.uniq)
        result = []
        if uniq < 50:
            result.append("low_map_rate")
        if unimap_gc > 45 or unimap_gc < 35:
            result.append("GC_unusual")
        if not len(result):
            result.append("-")
        self.content = ";".join(result)

    def count_lessdata(self):
        print "确定lessdata"
        if self.unimap_rds >= self.min_uniq_reads:
            self.lessdata = 'NO'
        else:
            self.lessdata = 'YES'


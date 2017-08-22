# -*- coding: UTF-8 -*-
import os
import glob
import fnmatch
import re
import shutil
#import time
#import numpy as np
#import rpy2.robjects as robjects
#from rpy2.robjects.packages import importr
from multiprocessing.dummy import Pool as ThreadPool
#from multiprocessing.pool import ThreadPool

class Trim(object):

    def __init__(self, nipt):
        self.indir = nipt.map_dir
        self.outdir = nipt.trim_dir
        self.gccorrect = nipt.gccorrect
        self.complexity = nipt.complexity
        self.merge = nipt.merge_block
        self.parallel = nipt.parallel
        self.sample = nipt.sample_deal
        self.bin_length = 20000
        dirlocate = os.path.dirname(os.path.abspath(__file__))
        self.rfile = '/'.join([dirlocate,'lib','R','pre_process.r'])

    def trim_pipeline(self):
        if not os.path.exists(self.indir):
            print ' '.join([self.indir, 'not exists'])
            exit(1)

        if not os.path.exists(self.outdir):
            #os.makedirs(self.outdir)
            pass
        else:
            old_staff = glob.glob(self.outdir)
            for staff in old_staff:
                shutil.rmtree(staff)

        rd_files = []
        gc_files = []
        pool = ThreadPool(self.parallel)
        for sample in self.sample:
            rd = glob.glob("/".join([self.indir, ''.join([sample,'*','20K.txt'])]))
            gc = glob.glob("/".join([self.indir, ''.join([sample,'*','20K.GC.txt'])]))
            if len(rd) == 1 and len(gc) == 1:
                rd_files.append(rd[0])
                gc_files.append(gc[0])
            else:
                raise Exception(' '.join([sample, 'rd or gc files not found or more than 1']))
        param_files = zip(rd_files, gc_files)
        pool.map(self.trim_single_sample, param_files)

    def trim_single_sample(self, args):
        rdfile, gcfile = args
        output_sampleprefix = '/'.join([self.outdir, rdfile.split('.')[0]])
        logfile = ''.join([output_sampleprefix, '.filter.log'])
        if os.path.exists(logfile): os.remove(logfile)
        # 第一种方案，直接调用R程序，为了减小难度
        command_list = [\
            '/usr/bin/Rscript', self.rfile,\
            '--gcfile', gcfile,\
            '--rcfile', rdfile,\
            '--outdir', self.outdir,\
            '--Merge', str(self.merge),\
        ]
        if self.gccorrect: command_list.append('--gccorrect')
        if self.complexity:
            command_list.append('--complexity')
            command_list.append(self.complexity.name)
        command = ' '.join(command_list)
        os.system(command)
        # 第二套方案把R程序嵌入到python中
        #utils = importr('utils')
        #base = importr('base')
        #gc = base.as_matrix(utils.read_table(file=gcfile,stringsAsFactors=False,row_names=1,sep='\t'))
        #rd = base.as_matrix(utils.read_table(file=rdfile,stringsAsFactors=False,row_names=1,sep='\t'))
        #print gc.dim
        #print rd.dim












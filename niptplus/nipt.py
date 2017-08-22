# -*- coding: UTF-8 -*-
import os
import glob
import re
import shutil
import niptplus.map_check
import niptplus.trim_data
import niptplus.pca
import niptplus.hmm

class NIPT(object):
    """一个NIPT类，用于每次的分析。实例记录了所有需要的信息，
    实例的方法对应不同的分析方式或者分析步骤。
    """

    def __init__(self, args):
        super(NIPT, self).__init__()
        if type(args) is dict :
            self.indir = args.get('gindir')
            self.outdir = args.get('goutdir')
            self.flowcell = args.get('flowcell')
            self.merge_block = args.get('merge')
            self.min_uniq_reads = args.get('min_uniq_reads')
            self.parallel = args.get('parallel')
            self._stage = args.get('stage')
            self.chrs = args.get('chr')
            self.gccorrect = args.get('gccorrect')
            self.complexity = args.get('complexity')
            self.map_dir = args.get('map_dir')
            self.trim_dir = args.get('trim_dir')
            self.pca_dir = args.get('pca_dir')
            self.hmm_dir = args.get('hmm_dir')
            self.hmm_zscore = args.get('hmm_zscore')
            self.hmm_length = args.get('hmm_length')
            self.hmm_rate = args.get('hmm_rate')
            self.hmm_zsd = args.get('hmm_zsd')
            self.hmm_exe = args.get('hmm_exe')
            self.hmm_disease = args.get('hmm_disease')
            self.step = args.get('step')
            self.rm_mean_pve = args.get('rm_mean_pve')
            self.rmPCs = args.get('rmPCs')
            self.redisseq = args.get('redisseq')
            self.pca_param_dir = args.get('pca_param_dir')
            self.sample_select = args.get('samples')
            self.__sample_candidate = set()
            self.__sample_deal = set()
        else:
            raise Exception('Input args should be dict')

    def get_stage(self):
        return self._stage

    def set_stage(self, stage):
        self._stage = stage
        if stage == 'check_map':
            rd_files = glob.glob('/'.join([self.indir, '*.20K.txt']))
            candidate = []
            for filename in rd_files:
                candidate.append(re.split('\.|_',re.split('/', filename)[-1])[0])
            self.__sample_candidate = set(candidate)
        elif stage == 'trim_data':
            rd_files = glob.glob('/'.join([self.map_dir, '*.20K.txt']))
            candidate = []
            for filename in rd_files:
                candidate.append(re.split('\.|_',re.split('/', filename)[-1])[0])
            self.__sample_candidate = set(candidate)
        elif stage == 'pca':
            files = glob.glob('/'.join([self.trim_dir, '*.filter.*.txt']))
            candidate = []
            for filename in files:
                candidate.append(re.split('\.|_', re.split('/', filename)[-1])[0])
            self.__sample_candidate = set(candidate)
        elif stage == 'hmm':
            matrix_files = glob.glob('/'.join([self.pca_dir, '*.normalized.zscore.matrix']))
            candidate = []
            with open(matrix_files[0], 'r') as f:
                f.readline()
                for line in f:
                    candidate.append(line.split("\t")[0])
            self.__sample_candidate = set(candidate)
        else:
            raise Exception("Stage must in check_map, trim_data, pca, hmm")

        if len(self.sample_select) == 0:
            self.__sample_deal = self.sample_candidate
        else:
            self.__sample_deal = self.sample_candidate.intersection(self.sample_select)
            diff_set = self.sample_select.difference(self.__sample_deal)
            if len(diff_set) > 0:
                raise Exception(" ".join([','.join(diff_set), 'not found']))

    stage = property(get_stage, set_stage)

    def set_sample_candidate(self, sample_set):
        if isinstance(sample_set, set):
            self.__sample_candidate = sample_set
        else:
            raise TypeError("Need set obj")
    def set_sample_deal(self, sample_set):
        if isinstance(sample_set, set):
            self.__sample_deal = sample_set
        else:
            raise TypeError("Need set obj")

    def get_sample_candidate(self):
        if len(self.__sample_candidate) > 0:
            return self.__sample_candidate
        if self.stage == 'check_map':
            rd_files = glob.glob('/'.join([self.indir, '*.20K.txt']))
            candidate = []
            for filename in rd_files:
                candidate.append(re.split('\.|_',re.split('/', filename)[-1])[0])
            self.__sample_candidate = set(candidate)
        elif self.stage == 'trim_data':
            rd_files = glob.glob('/'.join([self.map_dir, '*.20K.txt']))
            candidate = []
            for filename in rd_files:
                candidate.append(re.split('\.|_',re.split('/', filename)[-1])[0])
            self.__sample_candidate = set(candidate)
        elif self.stage == 'pca':
            files = glob.glob('/'.join([self.trim_dir, '*.filter.*.txt']))
            candidate = []
            for filename in files:
                candidate.append(re.split('\.|_', re.split('/', filename)[-1])[0])
            self.__sample_candidate = set(candidate)
        elif self.stage == 'hmm':
            matrix_files = glob.glob('/'.join([self.pca_dir, '*.normalized.zscore.matrix']))
            candidate = []
            with open(matrix_files[0], 'r') as f:
                f.readline()
                for line in f:
                    candidate.append(line.split("\t")[0])
            self.__sample_candidate = set(candidate)
        else:
            raise Exception("Stage must in check_map, trim_data, pca, hmm")

        if len(self.__sample_candidate) == 0:
            raise Exception('Sample candidate is empty')
        else:
            return self.__sample_candidate

    def get_sample_deal(self):
        if len(self.__sample_deal) > 0:
            return self.__sample_deal
        if len(self.sample_select) == 0:
            self.__sample_deal = self.sample_candidate
        else:
            self.__sample_deal = self.sample_candidate.intersection(self.sample_select)
            diff_set = self.sample_select.difference(self.__sample_deal)
            if len(diff_set) > 0:
                raise Exception(" ".join([','.join(diff_set), 'not found']))
        if len(self.__sample_deal) == 0:
            raise Exception('Sample needs to be deal is empty')
        return self.__sample_deal

    sample_candidate = property(get_sample_candidate, set_sample_candidate)
    sample_deal = property(get_sample_deal, set_sample_deal)

    def check_map(self):
        """检查对比结果
        """
        print "NIPT类进行比对检查"
        self.stage = 'check_map'
        check_obj = niptplus.map_check.CheckData(self)
        check_obj.check_map()
        print "check_map is done"

    def trim_data(self):
        print "NIPT类进行数据剪裁"
        self.stage = 'trim_data'
        trim_obj = niptplus.trim_data.Trim(self)
        trim_obj.trim_pipeline()
        print "trim_data is done"

    def pca_analyse(self):
        print "NIPT类进行主成分分析"
        #print vars(self)
        self.stage = 'pca'
        pca_obj = niptplus.pca.PCA(self)
        pca_obj.make_rd_matrix()
        pca_obj.pca()
        print "pca analyse is done"

    def hmm(self):
        print "NIPT类进行隐马分析"
        print vars(self)
        self.stage = 'hmm'
        hmm_obj = niptplus.hmm.HMM(self)
        print hmm_obj;
        hmm_obj.hmm_analyse();

    def pipeline(self):
        print "NIPT类的pipeline"
        if self.step == 'check_map':
            self.check_map()
            self.trim_data()
            self.pca_analyse()
            self.hmm()
        elif  self.step == 'trim_data':
            self.trim_data()
            self.pca_analyse()
            self.hmm()
        elif self.step == 'pca':
            self.pca_analyse()
            self.hmm()
        elif self.step == 'hmm':
            self.hmm()
        else:
            print "后面２个步骤还没有确定"

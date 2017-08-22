# -*- coding: UTF-8 -*-
import sys
import os
import niptplus
import niptplus.nipt
import argparse
import pprint
from inspect import getmembers

class SampleAction(argparse.Action):
    def __init__(self, option_strings, dest, nargs=None, **kwargs):
        if nargs is not None:
            raise ValueError("nargs not allowed")
        super(SampleAction, self).__init__(option_strings, dest, **kwargs)

    def __call__(self, parser, namespace, values, option_string=None):
        samples = set(values.split(','))
        setattr(namespace, self.dest, samples)

class ComplexityFileAction(argparse.Action):
    def __init__(self, option_strings, dest, nargs=None, **kwargs):
        super(ComplexityFileAction, self).__init__(option_strings, dest, nargs, **kwargs)
    def __call__(self, parser, namespace, values, option_string=None):
        if len(values) == 0:
            filename = '/'.join([os.path.dirname(os.path.abspath(__file__)),'data','refs','trim_data','20K_bins_complex.txt'])
            setattr(namespace, self.dest, open(filename, 'r'))
        else:
            setattr(namespace, self.dest, open(values[0], 'r'))

class ChrAction(argparse.Action):
    def __init__(self, option_strings, dest, nargs=None, **kwargs):
        super(ChrAction, self).__init__(option_strings, dest, nargs, **kwargs)
    def __call__(self, parser, namespace, values, option_string=None):
        if len(values) == 0:
            setattr(namespace, self.dest, set(xrange(1,23,1)))
        else:
            chr_str = values[0]
            chr_list = chr_str.split(',')
            chr_set = set()
            for chr_num in chr_list:
                chr_ = chr_num.split('.')
                if len(chr_) > 1:
                    chr_set = chr_set.union(set(xrange(int(chr_[0]),int(chr_[-1])+1,1 if chr_[1]=='' else int(chr_[1]))))
                else:
                    chr_set = chr_set.union(set([int(chr_num)]))
            setattr(namespace, self.dest, chr_set)

class PCAdirAction(argparse.Action):
    def __init__(self, option_strings, dest, nargs=None, **kwargs):
        super(PCAdirAction, self).__init__(option_strings, dest, nargs, **kwargs)
    def __call__(self, parser, namespace, values, option_string=None):
        if len(values) == 0:
            filedir = '/'.join([os.path.dirname(os.path.abspath(__file__)),'data','refs','pca_data/refs',])
            setattr(namespace, self.dest, filedir)
        else:
            setattr(namespace, self.dest, values[0])


class PCAredisseqAction(argparse.Action):
    def __init__(self, option_strings, dest, nargs=None, **kwargs):
        super(PCAredisseqAction, self).__init__(option_strings, dest, nargs, **kwargs)
    def __call__(self, parser, namespace, values, option_string=None):
        if len(values) == 0:
            filedir = '/'.join([os.path.dirname(os.path.abspath(__file__)),'data','refs','pca_data/20K.CHRs.block',])
            setattr(namespace, self.dest, filedir)
        else:
            setattr(namespace, self.dest, values[0])


def niptplus_commands():
    parser = argparse.ArgumentParser()
    general_option = parser.add_argument_group("General Options")
    general_option.add_argument('--gindir', help='Input base dir', default='/seq_dir/prenatal')
    general_option.add_argument('--goutdir', help='Output base dir', default='.')
    general_option.add_argument('--flowcell', '-f', help='Flowcell to analyse', required=True)
    general_option.add_argument('--samples', '-s', help='Sample to analyse', default=set(), action=SampleAction)
    general_option.add_argument('--merge', '-m', help='Merge block to use', default=5, type=int)
    general_option.add_argument('--min_uniq_reads', '-min', help='Min reads for pass', type=int, default=11000000)
    general_option.add_argument('--parallel', '-par', help='Number of parallel', type=int, default=10)
    general_option.add_argument('--chr', '-chr', help='Chr to analyse', action=ChrAction, default=set(xrange(1,23,1)), nargs=1)
    subparsers = parser.add_subparsers(title='Commands', description='valid subcommands', help='sub-command help')

    parser_echo = subparsers.add_parser('echo', help="Print str which you use")
    parser_echo.add_argument("echo", help="echo the string you use here")
    parser_echo.set_defaults(func=echo)

    parser_fun = subparsers.add_parser('fun', help="Try fun function")
    parser_fun.set_defaults(func=fun)

    parser_step1 = subparsers.add_parser('check_map', help='Check mapping result')
    parser_step1.set_defaults(func=check_map)
    parser_step1.add_argument("--indir", metavar="GINDIR/FLOWCELL", help="Input dir which stores mapping result")
    parser_step1.add_argument("--outdir", metavar="GOUTDIR/map_data", help="Output dir for checking result")

    parser_step2 = subparsers.add_parser('trim_data', help="Trim data to match rules")
    parser_step2.set_defaults(func=trim_data)
    parser_step2.add_argument("--indir", metavar="GOUTDIR/map_data", help="Input dir for checked mapping result")
    parser_step2.add_argument("--outdir", metavar="GOUTDIR/trim_data", help="Output dir for triming result")
    parser_step2.add_argument('--complexity', help='Filter bins with complexity rate of each bins', action=ComplexityFileAction, nargs='*')
    parser_step2.add_argument('--gccorrect', help='GC correction', action='store_true')

    parser_step3 = subparsers.add_parser('pca', help="PCA analysis for trimed data")
    parser_step3.set_defaults(func=pca_analyse)
    parser_step3.add_argument("--indir", metavar="GOUTDIR/trim_data", help="Indir for trimed data")
    parser_step3.add_argument("--outdir", metavar="GOUTDIR/pca_data", help="outdir for pca result")
    parser_step3.group = parser_step3.add_mutually_exclusive_group()
    parser_step3.group.add_argument("--rm_mean_pve", '-rm_m_p', const=0.95, action='store', help="Unknown", nargs="?")
    parser_step3.group.add_argument('--rmPCs','-rmpcs', help='the move out PCs number', action='store', const=10, nargs="?")
    parser_step3.add_argument('--redisseq', action=PCAredisseqAction, nargs='*')
    parser_step3.add_argument('--param_dir',dest='pca_param_dir', action=PCAdirAction, nargs='*')


    parser_step4 = subparsers.add_parser('hmm', help="HMM for pca data")
    parser_step4.set_defaults(func=hmm)
    parser_step4.add_argument("--indir", metavar="GOUTDIR/pca_data", help="Indir for hmm")
    parser_step4.add_argument("--outdir", metavar="GOUTDIR/hmm_data", help="Outdir for hmm")
    parser_step4.add_argument("--hmm_zscore", '-h_zs', default=2, help="cutoff of zscore in observations of HMM, default 2", type=float)
    parser_step4.add_argument("--hmm_length", '-h_len', default=6, help="the mean number of bins in a cnv, default 6", type=float)
    parser_step4.add_argument("--hmm_rate", '-h_rate', default='1e-5', help="the chromosome wide rate, default 1e-5")
    parser_step4.add_argument("--hmm_zsd", '-h_zsd', default=1, help="the Standard deviation of del/dup z-score distribution, default 1", type=float)
    parser_step4.add_argument('--hmm_exe', '-h_exe',
            default='/'.join([os.path.dirname(os.path.abspath(__file__)),'lib','xhmm']),
            help="The path for hmm exe")
    parser_step4.add_argument('--hmm_disease', '-h_disease', help="Ref for 7 disease",
            default='/'.join([os.path.dirname(os.path.abspath(__file__)),'data','refs','hmm/disease.txt',]))

    parser_pipeline = subparsers.add_parser('pipeline', help="Pipeline for analyse")
    parser_pipeline.set_defaults(func=pipeline)
    parser_pipeline.add_argument("--step", help="Step for start", default="check_map", choices=["check_map", "trim_data", "pca", "hmm", "visual", 'summary'])
    parser_pipeline.add_argument("--map_dir", help="\tResult dir of mapping check", default='map_data')
    parser_pipeline.add_argument("--trim_dir", help="\tResult dir of triming data", default="trim_data")
    parser_pipeline.add_argument("--pca_dir", help="\tResult dir of pca", default="pca_data")
    parser_pipeline.add_argument("--hmm_dir", help="\tResult dir of hmm", default="hmm_data")

    parser_pipeline_trim_stage = parser_pipeline.add_argument_group("Trim stage options")
    parser_pipeline_trim_stage.add_argument('--complexity', help='Filter bins with complexity rate of each bins', action=ComplexityFileAction, nargs='*')
    parser_pipeline_trim_stage.add_argument('--gccorrect', help='GC correction', action='store_true')

    parser_pipeline_pca_stage = parser_pipeline.add_argument_group("PCA stage options")
    parser_pipeline_pca_stage.group = parser_pipeline_pca_stage.add_mutually_exclusive_group()
    parser_pipeline_pca_stage.group.add_argument("--rm_mean_pve", '-rm_m_p', const=0.95, action='store', help="Unknown", nargs="?")
    parser_pipeline_pca_stage.group.add_argument('--rmPCs','-rmpcs', help='the move out PCs number', action='store', const=10, nargs="?")
    parser_pipeline_pca_stage.add_argument('--redisseq', action=PCAredisseqAction, nargs='*')
    parser_pipeline_pca_stage.add_argument('--param_dir',dest='pca_param_dir', action=PCAdirAction, nargs='*')

    parser_pipeline_hmm_stage = parser_pipeline.add_argument_group("HMM stage options")
    parser_pipeline_hmm_stage.add_argument("--hmm_zscore", '-h_zs', default=2, help="cutoff of zscore in observations of HMM, default 2", type=float)
    parser_pipeline_hmm_stage.add_argument("--hmm_length", '-h_len', default=6, help="the mean number of bins in a cnv, default 6", type=float)
    parser_pipeline_hmm_stage.add_argument("--hmm_rate", '-h_rate', default='1e-5', help="the chromosome wide rate, default 1e-5")
    parser_pipeline_hmm_stage.add_argument("--hmm_zsd", '-h_zsd', default=1, help="the Standard deviation of del/dup z-score distribution, default 1", type=float)
    parser_pipeline_hmm_stage.add_argument('--hmm_exe', '-h_exe',
        default='/'.join([os.path.dirname(os.path.abspath(__file__)),'lib','xhmm']),
        help="The path for hmm exe")
    parser_pipeline_hmm_stage.add_argument('--hmm_disease', '-h_disease', help="Ref for 7 disease",
        default='/'.join([os.path.dirname(os.path.abspath(__file__)),'data','refs','hmm/disease.txt',]))


    if len(sys.argv[1:])==0:
        parser.print_help()
        return

    args = parser.parse_args()
    args.goutdir = '/'.join([args.goutdir, args.flowcell])
    args.gindir = "/".join([args.gindir, args.flowcell])
    args.func(args)

def set_indir_and_outdir(args, basein, indir, baseout, outdir):
    if not args.indir:
        args.indir = '/'.join([basein, indir])
    if not args.outdir:
        args.outdir = '/'.join([baseout, outdir])

def echo(args):
    print args.echo

def fun(args):
    print "This is niptplus commands fun function"
    print niptplus.joke()
    print "-"*30
    print vars(args)
    print "-"*30
    pprint.pprint(getmembers(args))

def check_map(args):
    #print "开发check map"
    set_indir_and_outdir(args, args.gindir, args.flowcell, args.goutdir, 'map_data')
    #print vars(args)
    args.__setattr__('map_dir', args.outdir)
    #args.__setattr__('stage', 'check_map')
    obj = niptplus.nipt.NIPT(args.__dict__)
    obj.check_map()
    #print vars(obj)
    #print "check map开发完成"

def trim_data(args):
    #print "开发trim data"
    base = args.goutdir
    set_indir_and_outdir(args, base, 'map_data', base, 'trim_data')
    args.__setattr__('map_dir', args.indir)
    args.__setattr__('trim_dir', args.outdir)
    #args.__setattr__('stage', 'trim_data')
    obj = niptplus.nipt.NIPT(args.__dict__)
    obj.trim_data()
    #print vars(obj)
    #print "trim_data开发完成"

def pca_analyse(args):
    #print "开发pca analyse"
    base = args.goutdir
    set_indir_and_outdir(args, base, 'trim_data', base, 'pca_data')
    args.__setattr__('trim_dir', args.indir)
    args.__setattr__('pca_dir', args.outdir)
    #args.__setattr__('stage', 'pca')
    obj = niptplus.nipt.NIPT(args.__dict__)
    obj.pca_analyse()
    #print vars(obj)
    #print "pca analyse开发完成"

def hmm(args):
    print "开发hmm"
    base = args.goutdir
    set_indir_and_outdir(args, base, 'pca_data', base, 'hmm_data')
    args.__setattr__('pca_dir', args.indir)
    args.__setattr__('hmm_dir', args.outdir)
    #args.__setattr__('stage', 'hmm')
    #print vars(args)
    obj = niptplus.nipt.NIPT(args.__dict__)
    obj.hmm()
    #print vars(obj)

def pipeline(args):
    print "开发pipeline"
    for attr in ["map_dir", "trim_dir", 'pca_dir', 'hmm_dir'] :
        args.__setattr__(attr, '/'.join([args.goutdir, args.__getattribute__(attr)]))
    #args.__setattr__('stage', args.get('step'))
    #print vars(args)
    obj = niptplus.nipt.NIPT(args.__dict__)
    obj.pipeline()
    #pprint.pprint(getmembers(args))
    #pprint.pprint(getmembers(obj))


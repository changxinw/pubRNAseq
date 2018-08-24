import os
import re
import sys

def check_design(dmat):
    if not os.path.exists(dmat):
        os.system("echo Design matrix not exist!")
        sys.exit(0)
        return 0
    with open(dmat) as design_mat:
        for i in design_mat.realines():
            if len(i.rstrip().split('\t')) != 2:
                os.system("echo Wrong design matrix!")
                sys.exit(0)
                return 0
            elif not re.match(r'^GSM\d+$', i.rstrip().split('\t')[0]):
                os.system("echo Wrong GSM ID provide!")
                sys.exit(0)
                return 0
            else:
                return 1

def opt_check(opt):
    """
    Check each option and the format of design matrix
    return true or false
    """
    if not opt.dmat:
        os.system("echo Please give a design matrix with GSM ID and experiment design!")
        sys.exit(0)
    elif not opt.output:
        os.system("echo Please give a directory to output result!")
        sys.exit(0)
    elif not os.path.isdir(opt.output):
        os.system("echo Output path not exist, pyRNAseq make it for you!")
        os.makedirs(opt.output)
    elif opt.salmon and not opt.species:
        os.system("echo Please define your species as hg38 or mm10 if you need mapping!")
        sys.exit(0)
    elif opt.difexpr and not opt.salmon:
        os.system("echo Please using salmon to align if your want to differential expression analysis!")
        sys.exit(0)
    else:
        return 1
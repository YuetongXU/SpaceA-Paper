#!/usr/bin/env python
# split_dpm3_dpm4_from_full_BC_fq

# encoding: utf-8
# Split fq with index
# conda activate myenv

import os
import argparse
import re
from Bio.SeqIO.QualityIO import FastqGeneralIterator
from gzip import open as gzopen

def parse_arguments():
    parser = argparse.ArgumentParser(description='split_dpm3_dpm4_from_full_BC_fq')
    parser.add_argument('--full_bc_fq',dest='full_bc_fq', required=True,help='full_BC_fq')
    parser.add_argument('--dpm3',dest='dpm3', required=True,help='dpm3')
    parser.add_argument('--dpm4',dest='dpm4', required=True,help='dpm4')
    return parser.parse_args()

def main(args):
    full_bc_fq=os.path.join(os.getcwd(),args.full_bc_fq)
    dpm3=os.path.join(os.getcwd(),args.dpm3)
    dpm4=os.path.join(os.getcwd(),args.dpm4)

    total_num=0
    dpm3_num=0
    dpm4_num=0
    others_num=0
    pattern = re.compile('\[([a-zA-Z0-9_\-]+)\]')
    
    with gzopen(full_bc_fq,'rt') as in_h,gzopen(dpm3,'wt') as dpm3_h,gzopen(dpm4,'wt') as dpm4_h:
        for title, seq, qual in FastqGeneralIterator(in_h):
            total_num+=1
            barcodes = pattern.findall(title)

            if 'dpm3' in barcodes:
                dpm3_h.write("@%s\n%s\n+\n%s\n" % (title, seq, qual))
                dpm3_num+=1
            elif 'dpm4' in barcodes:
                dpm4_h.write("@%s\n%s\n+\n%s\n" % (title, seq, qual))
                dpm4_num+=1
            else:
                others_num+=1

        print('total_num: ',total_num)
        print('dpm3_num: ',dpm3_num)
        print('dpm4_num: ',dpm4_num)
        print('others_num: ',others_num)

        print('dpm3_num/total_num(%): ',round(dpm3_num/total_num*100,2))
        print('dpm4_num/total_num(%): ',round(dpm4_num/total_num*100,2))
        print('others_num/total_num(%): ',round(others_num/total_num*100,2))

if __name__ == "__main__":
    args =parse_arguments()
    main(args)




import math
import sys
import os
from collections import defaultdict


def nCr(n, r):
    'cal contact number'
    f = math.factorial
    return f(n) // f(r) // f(n - r)


def get_chr_read_dict(C_sprite_read: list) -> dict:
    chr_read_dict = defaultdict(list)
    for read in C_sprite_read:
        chr = read.split(']_')[1].split(':')[0]
        chr_read_dict[chr].append(read)
    return chr_read_dict


def get_chr_read_num_dict(chr_read_dict: dict) -> dict:
    return {chr: len(chr_read_dict[chr]) for chr in chr_read_dict.keys()}


def get_cis_ct_num(chr_read_num_dict: dict) -> int:
    return sum((nCr(x, 2) for x in chr_read_num_dict.values() if x >= 2))


def get_chr_read_num_ratio_dict(chr_read_num_dict: dict, total_read_num: int) -> dict:
    return {chr: round(x / total_read_num * 100, 2) for chr, x in chr_read_num_dict.items()}


def spot_qc(spot_id: str, sample_id: str, spot_f_dir: str, min_C_size: int, max_C_size: int,
            min_chr_read_num_ratio: float, min_C_cis_ct_num_ratio: float):
    
    spot_infor = {}
    Cs_infor = {}

    spot_raw_read_num = 0
    spot_fcsize_read_num = 0
    spot_raw_C_num = 0
    spot_fcsize_C_num = 0
    spot_raw_ct_num = 0
    spot_fcsize_ct_num = 0
    spot_fcsizefchrrnp_ct_num = 0
    spot_fcsize_cis_ct_num = 0
    spot_fcsizefchrrnp_cis_ct_num = 0
    spot_fcsizefchrrnp_read_num = 0

    spot_sprite_f_path = os.path.join(spot_f_dir, f'clusters_{sample_id}_single_filtered/{spot_id}')
    with open(spot_sprite_f_path, 'r') as f:
        # cluster
        for line in f:
            spot_raw_C_num += 1
            reads = line.rstrip().split()[1:]
            C_id = line.rstrip().split()[0]

            unique_read = set(reads)
            spot_raw_read_num += len(unique_read)
            C_size = len(unique_read)

            if C_size >= 2:
                spot_raw_ct_num += nCr(C_size, 2)
            else:
                continue

            # filter cluster of cluster size
            if min_C_size <= C_size <= max_C_size:
                C_infor = {}
                C_infor['fcsize_C_id'] = C_id
                C_infor['fcsize_size'] = C_size

                spot_fcsize_C_num += 1
                spot_fcsize_read_num += C_size
                C_fcsize_ct_num = nCr(C_size, 2)
                C_infor['fcsize_ct_num'] = C_fcsize_ct_num

                spot_fcsize_ct_num += C_fcsize_ct_num

                C_fcsize_chr_read_dict = get_chr_read_dict(unique_read)
                C_infor['fcsize_chr_read'] = C_fcsize_chr_read_dict

                C_fcsize_chr_read_num_dict = get_chr_read_num_dict(C_fcsize_chr_read_dict)
                C_fcsize_cis_ct_num = get_cis_ct_num(C_fcsize_chr_read_num_dict)
                C_infor['fcsize_cis_ct_num'] = C_fcsize_cis_ct_num

                C_infor['fcsize_cis_ct_ratio'] = C_fcsize_cis_ct_num / C_fcsize_ct_num if C_fcsize_ct_num != 0 else 0

                spot_fcsize_cis_ct_num += C_fcsize_cis_ct_num

                # filter cluster reads of cluster chr reads num pct low
                C_fcsize_chr_read_num_ratio_dict = {
                    chr: C_fcsize_chr_read_num_dict[chr] / C_size for chr in C_fcsize_chr_read_num_dict.keys()}
                C_infor['fcsize_chr_read_num_ratio'] = C_fcsize_chr_read_num_ratio_dict

                C_fcsizefchrrnp_chrs = [
                    chr for chr in C_fcsize_chr_read_num_ratio_dict.keys() if
                    C_fcsize_chr_read_num_ratio_dict[chr] > min_chr_read_num_ratio]
                C_infor['fcsizefchrrnp_chrs'] = C_fcsizefchrrnp_chrs
                C_fcsizefchrrnp_chr_read_dict = {
                    chr: C_fcsize_chr_read_dict[chr] for chr in C_fcsizefchrrnp_chrs}
                C_infor['fcsizefchrrnp_chr_read'] = C_fcsizefchrrnp_chr_read_dict

                C_fcsizefchrrnp_chr_read_num_dict = get_chr_read_num_dict(C_fcsizefchrrnp_chr_read_dict)
                C_fcsizefchrrnp_read_num = sum(C_fcsizefchrrnp_chr_read_num_dict.values())
                C_infor['fcsizefchrrnp_read_num'] = C_fcsizefchrrnp_read_num
                spot_fcsizefchrrnp_read_num += C_fcsizefchrrnp_read_num

                C_fcsizefchrrnp_ct_num = nCr(C_fcsizefchrrnp_read_num, 2) if C_fcsizefchrrnp_read_num >= 2 else 0
                C_infor['fcsizefchrrnp_ct_num'] = C_fcsizefchrrnp_ct_num
                spot_fcsizefchrrnp_ct_num += C_fcsizefchrrnp_ct_num

                C_fcsizefchrrnp_cis_ct_num = get_cis_ct_num(C_fcsizefchrrnp_chr_read_num_dict)
                C_infor['fcsizefchrrnp_cis_ct_num'] = C_fcsizefchrrnp_cis_ct_num
                spot_fcsizefchrrnp_cis_ct_num += C_fcsizefchrrnp_cis_ct_num

                C_fcsizefchrrnp_cis_ct_num_ratio = C_fcsizefchrrnp_cis_ct_num / C_fcsizefchrrnp_ct_num if C_fcsizefchrrnp_ct_num != 0 else 0
                C_infor['fcsizefchrrnp_cis_ct_num_ratio'] = C_fcsizefchrrnp_cis_ct_num_ratio

                # cluster cis contact num pct
                C_cis_ct_num_ratio_low = True if C_fcsizefchrrnp_cis_ct_num_ratio < min_C_cis_ct_num_ratio else False
                C_infor['C_cis_ct_num_ratio_low'] = C_cis_ct_num_ratio_low
                Cs_infor[C_id] = C_infor

                # 释放不再使用的中间变量
                del C_fcsize_chr_read_dict, C_fcsize_chr_read_num_dict, C_fcsizefchrrnp_chr_read_dict, C_fcsizefchrrnp_chr_read_num_dict

        #  filter cluster of cluster cis contact num pct low
        Cs_infor_fCsofcis = {C_id: C_infor for C_id, C_infor in Cs_infor.items() if
                             C_infor['C_cis_ct_num_ratio_low'] == False}

        spot_fcsizefchrrnpfCsofcis_C_num = len(Cs_infor_fCsofcis.keys())
        spot_fcsizefchrrnpfCsofcis_read_num = sum(
            [Cs_infor_fCsofcis[C_id]['fcsizefchrrnp_read_num'] for C_id in Cs_infor_fCsofcis.keys()])
        spot_fcsizefchrrnpfCsofcis_ct_num = sum(
            [Cs_infor_fCsofcis[C_id]['fcsizefchrrnp_ct_num'] for C_id in Cs_infor_fCsofcis.keys()])
        spot_fcsizefchrrnpfCsofcis_cis_ct_num = sum(
            [Cs_infor_fCsofcis[C_id]['fcsizefchrrnp_cis_ct_num'] for C_id in Cs_infor_fCsofcis.keys()])
        spot_fcsizefchrrnpfCsofcis_cis_ct_num_ratio = spot_fcsizefchrrnpfCsofcis_cis_ct_num / spot_fcsizefchrrnpfCsofcis_ct_num if spot_fcsizefchrrnpfCsofcis_ct_num != 0 else 0

        spot_infor['spot_raw_read_num'] = spot_raw_read_num
        spot_infor['spot_raw_C_num'] = spot_raw_C_num
        spot_infor['spot_raw_ct_num'] = spot_raw_ct_num
        spot_infor['spot_fcsize_read_num'] = spot_fcsize_read_num
        spot_infor['spot_fcsize_C_num'] = spot_fcsize_C_num
        spot_infor['spot_fcsize_ct_num'] = spot_fcsize_ct_num
        spot_infor['spot_fcsize_cis_ct_num'] = spot_fcsize_cis_ct_num
        spot_infor['spot_fcsize_cis_ct_num_ratio'] = spot_fcsize_cis_ct_num / spot_fcsize_ct_num if spot_fcsize_ct_num != 0 else 0
        spot_infor['spot_fcsizefchrrnp_read_num'] = spot_fcsizefchrrnp_read_num
        spot_infor['spot_fcsizefchrrnp_ct_num'] = spot_fcsizefchrrnp_ct_num
        spot_infor['spot_fcsizefchrrnp_cis_ct_num'] = spot_fcsizefchrrnp_cis_ct_num
        spot_infor['spot_fcsizefchrrnp_cis_ct_num_ratio'] = spot_fcsizefchrrnp_cis_ct_num / spot_fcsizefchrrnp_ct_num if spot_fcsizefchrrnp_ct_num != 0 else 0
        spot_infor['spot_fcsizefchrrnpfCsofcis_C_num'] = spot_fcsizefchrrnpfCsofcis_C_num
        spot_infor['spot_fcsizefchrrnpfCsofcis_read_num'] = spot_fcsizefchrrnpfCsofcis_read_num
        spot_infor['spot_fcsizefchrrnpfCsofcis_ct_num'] = spot_fcsizefchrrnpfCsofcis_ct_num
        spot_infor['spot_fcsizefchrrnpfCsofcis_cis_ct_num'] = spot_fcsizefchrrnpfCsofcis_cis_ct_num
        spot_infor['spot_fcsizefchrrnpfCsofcis_cis_ct_num_ratio'] = spot_fcsizefchrrnpfCsofcis_cis_ct_num_ratio

        # 释放不再使用的中间变量
        del Cs_infor
        # print(f'{spot_id}',end=' ',sep='')
    return spot_infor, Cs_infor_fCsofcis, Cs_infor_fCsofcis
def spot_filtered_Cs_infor_to_sprite(Cs_infor,out_sprite_path):
    with open(out_sprite_path, 'w') as f:
        for cluster in Cs_infor.keys():
            chr_read_dict=Cs_infor[cluster]['fcsizefchrrnp_chr_read']
            C_reads=[]
            for chr in chr_read_dict.keys():
                chr_reads=chr_read_dict[chr]
                C_reads+=chr_reads
            f.write('\t'.join([cluster]+C_reads)+'\n')

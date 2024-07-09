#!/usr/bin/env python
import sys
import numpy as np
import pandas as pd
import re
import itertools
import math
import multiprocessing as mp
from threading import Thread
import Bio
from Bio import SeqIO
import os
import matplotlib.pyplot as plt
from plotnine import *
import seaborn as sns
from sklearn.linear_model import LinearRegression
#import pymol
#import umap
#import umap.plot


script_path = os.path.abspath(os.path.join(os.path.dirname(os.path.realpath(__file__)), '../bbmap/current/'))
## check if minimap2 is installed
if os.system('minimap2 --version') != 0:
    minimap_path = os.path.abspath(os.path.join(os.path.dirname(os.path.realpath(__file__)), '../minimap2-2.26_x64-linux/minimap2 '))
else:
    minimap_path = 'minimap2 '

codon_chart = {'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L', 'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L', 'ATT': 'I',
         'ATC': 'I', 'ATA': 'I', 'ATG': 'M', 'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V', 'TCT': 'S', 'TCC': 'S',
         'TCA': 'S', 'TCG': 'S', 'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P', 'ACT': 'T', 'ACC': 'T', 'ACA': 'T',
         'ACG': 'T', 'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A', 'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
         'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q', 'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K', 'GAT': 'D',
         'GAC': 'D', 'GAA': 'E', 'GAG': 'E', 'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W', 'CGT': 'R', 'CGC': 'R',
         'CGA': 'R', 'CGG': 'R', 'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R', 'GGT': 'G', 'GGC': 'G', 'GGA': 'G',
         'GGG': 'G'}


def merge_pairs(in1, in2, out, bbmerge_location=f'java -cp {script_path} jgi.BBMerge'):
    bbmerge_command = bbmerge_location
    merge_out = os.path.splitext(in1)[0]+'_merge'+os.path.splitext(in1)[1]
    unmerge_out1 = os.path.splitext(in1)[0]+'_unmerge'+os.path.splitext(in1)[1]
    unmerge_out2 = os.path.splitext(in2)[0]+'_unmerge'+os.path.splitext(in2)[1]
    for file in [out, merge_out, unmerge_out1, unmerge_out2]:
        if os.path.exists(file):
            os.remove(file)
    bbmerge_command += f' in1={in1} in2={in2} out={merge_out} outu1={unmerge_out1} outu2={unmerge_out2}'
    print(f'{bbmerge_command}')
    ret = os.system(bbmerge_command)
    assert ret == 0
    filenames = [merge_out, unmerge_out1, unmerge_out2]
    with open(out, 'w') as outfile:
        for fname in filenames:
            with open(fname) as infile:
                for line in infile:
                    outfile.write(line)
    for tempfile in [merge_out, unmerge_out1, unmerge_out2]:
        os.remove(tempfile)


def correct_pairs(in1, in2, bbmerge_location=f'java -cp {script_path} jgi.BBMerge'):
    bbmerge_command = bbmerge_location
    if 'fastq.gz' in in1:
        out1 = in1.split('.fastq.gz')[0]+'_corrected.fastq.gz'
        out2 = in2.split('.fastq.gz')[0] + '_corrected.fastq.gz'
    else:
        out1 = os.path.splitext(in1)[0]+'_corrected'+os.path.splitext(in1)[1]
        out2 = os.path.splitext(in2)[0]+'_corrected'+os.path.splitext(in2)[1]
    bbmerge_command += f' in1={in1} in2={in2} out1={out1} ecco mix'
    print(bbmerge_command)
    ret = os.system(bbmerge_command)
    assert ret == 0
    return out1


def chunkit(align_files, file_size, chunk_size, other_args=()):
    # define functions for finding breakpoints
    def is_start_of_line(position):
        if position == 0:
            return True
        # Check whether the previous character is end of line
        f.seek(position - 1)
        return f.read(1) == '\n'

    def get_next_line_position(position):
        # Read the current line till the end
        f.seek(position)
        f.readline()
        # Return a position after reading the line
        return f.tell()

    def is_paired(position):
        # Read the current line till the end
        f.seek(position)
        name = f.readline().split('\t')[0]
        # Return a position after reading the line
        return next(f).split('\t')[0] == name

    # Arguments for each chunk
    chunk_args = []
    with open(align_files, 'r') as f:
        # Iterate over all chunks and construct arguments for `process_chunk`
        chunk_start = 0
        while chunk_start < file_size:
            chunk_end = min(file_size, chunk_start + chunk_size)
            # Make sure the chunk ends at the beginning of the next line
            while not is_start_of_line(chunk_end):
                chunk_end -= 1
            # If the next line is paired, move the chunk end to the next line
            f.seek(chunk_end)
            if (chunk_end + len(f.readline())) < file_size:
                while not is_paired(chunk_end):
                    f.seek(chunk_end)
                    chunk_end += len(f.readline())
                while not is_start_of_line(chunk_end):
                    chunk_end -= 1
            # Handle the case when a line is too long to fit the chunk size
            if chunk_start == chunk_end:
                chunk_end = get_next_line_position(chunk_end)
            # Save `process_chunk` arguments
            args = (align_files, chunk_start, chunk_end) + other_args
            chunk_args.append(args)
            # Move to the next chunk
            chunk_start = chunk_end
    return chunk_args


def david_call_parallel_variants(sam_file, wt, outfile, app):
    cpu_count = min(app.parallel, mp.cpu_count())
    print('Running parallel analysis with', cpu_count, 'cores')
    file_size = os.path.getsize(sam_file)
    chunk_size = file_size // cpu_count
    chunk_args = chunkit(sam_file, file_size, chunk_size, other_args=(wt, outfile, app))
    # parallel
    with mp.Pool(cpu_count) as pool:
        answers = pool.starmap(david_call_variants, chunk_args)
    mutation_dict = {}
    variants = {'WT': 0}
    wt_count = [0]*int(len(wt)/3)
    for answer in answers:
        mutation_dict.update(answer[0])
        variants.update(answer[1])
        wt_count = [sum(x) for x in zip(wt_count, answer[2])]
    write_output(mutation_dict, variants, wt_count, wt, outfile)


def parse_cigar(r3, wt, quality_nt, tmpcounts, read_length, tmp_wt_count, indel_size, bowtie):
    # r3 is the cigar string split into a list
    # 0 is the read name
    # 1 is the flag (0 is forward, 16 is reverse)
    # 2 is the reference name
    # 3 is the reference position
    # 4 is the mapping quality
    # 5 is the cigar string
    # 6 is the mate reference name
    # 7 is the mate reference position
    # 8 is the template length
    # 9 is the read sequence
    # 10 is the read quality
    # 11 is the optional fields
    # 12 is the optional fields
    # set index
    ref = int(r3[3]) - 1  # subtract 1 for python numbering
    read_length[0] = min(read_length[0], ref)
    read = 0
    # cycle through each codon until end of read or end of reference sequencee
    fragment_cigar = list(filter(None, re.split(r'(\d+)', r3[5])))
    for cigar in [fragment_cigar[ii:ii + 2] for ii in range(0, len(fragment_cigar), 2)]:
        if 'X' == cigar[1] or ('M' == cigar[1] and bowtie == 1):
            # this might be a substitution
            for i_sub in range(int(cigar[0])):  # not sure how this range works
                codon_diff = (ref + i_sub) % 3  # find the correct frame
                i_codon = r3[9][read + i_sub - codon_diff:read + i_sub - codon_diff + 3].upper()
                wt_codon = str(wt.seq[ref + i_sub - codon_diff: ref + i_sub - codon_diff + 3])
                # record if codons are different
                if i_codon != wt_codon:
                    # record codon if you can see the entire codon
                    if read + i_sub - codon_diff + 3 < len(
                            r3[9]) and read + i_sub - codon_diff >= 0 and 'N' not in i_codon and '-' not in i_codon:
                        # only record if confident about entire codon
                        quality_codon = r3[10][read + i_sub - codon_diff:read + i_sub - codon_diff + 3]
                        if all([ord(x) - 33 > quality_nt for x in quality_codon]):
                            if codon_chart[i_codon] == codon_chart[wt_codon] and not bowtie:
                                tmpcounts.add(str(int((ref + i_sub - codon_diff) / 3 + 1)) + '_' + i_codon + '*')
                            else:
                                tmpcounts.add(str(int((ref + i_sub - codon_diff) / 3 + 1)) + '_' + i_codon)
            ref += int(cigar[0])
            read += int(cigar[0])
        elif '=' == cigar[1]:
            # this is a perfect match
            codon_start = ref + ref % 3
            codon_end = ref + int(cigar[0]) - (int(cigar[0]) % 3)
            if codon_end - codon_start >= 3:  # if the perfect match spans an entire codon then count it
                # make sure this doesn't exceed the length of the gene using min
                wt_range = range(int(codon_start / 3 + 1), min(int((codon_end - 3) / 3 + 1) + 1, int(len(wt) / 3)),
                                 1)  # need to add 1 for range
                for wt_i in wt_range:
                    tmp_wt_count[wt_i] = 1
            ref += int(cigar[0])
            read += int(cigar[0])
        elif 'D' == cigar[1]:
            # this is a deletion
            tmpcounts.add(str(int(ref / 3 + 1)) + '_' + 'del' + cigar[0])
            # insert '-' for gap
            r3[9] = r3[9][:read] + '-' * int(cigar[0]) + r3[9][read:]
            r3[10] = r3[10][:read] + '-' * int(cigar[0]) + r3[9][read:]
            read += int(cigar[0])
            ref += int(cigar[0])
            indel_size += int(cigar[0])
        elif 'I' == cigar[1]:
            # this is an insertion
            insertion = r3[9][read:read + int(cigar[0])].upper()
            tmpcounts.add(str(int(ref / 3 + 1)) + '_' + 'ins' + insertion)
            read += int(cigar[0])
            indel_size += int(cigar[0])
        elif 'S' == cigar[1]:
            # soft clipping
            read += int(cigar[0])
        elif 'M' == cigar[1] and not bowtie:
            # bad read (N) - ignore
            ref += int(cigar[0])
            read += int(cigar[0])
        elif 'H' == cigar[1]:
            # hard clipping
            pass
        elif 'P' == cigar[1]:
            # padding
            pass
        elif 'N' == cigar[1]:
            # skipped region from the reference
            ref += int(cigar[0])
        else:
            raise TypeError('Cigar format not found')
    read_length[1] = max(read_length[1], ref)
    return [tmpcounts, read_length, tmp_wt_count, indel_size]


def david_call_variants(sam_file, start, end, wt, outfile, app):
    quality = int(app.minq)
    quality_nt = int(app.minb)
    # Analyze Codons
    loop_read = open(sam_file, 'r')
    loop_read.seek(start)
    mutation_dict = {}
    variants = {'WT': 0}
    variant_tracking = {'WT': {'read_names': []}}
    wt_count = [0]*int(len(wt)/3)
    file = open(outfile + '_' + str(start) + '.csv', 'w')
    wt_file = open(outfile + '_' + str(start) + '_wt.csv', 'w')
    variant_check = app.variant
    correlation_check = app.correlation
    bowtie = 0
    while True:  # loop through each read from ngs
        try:
            line = next(loop_read)
            read_position = start
            tmp_r1 = line.strip().split('\t')
            start += len(line)
            if start > end:
                break
        except:
            break
        if '@' in tmp_r1[0]:
            if any(1 for x in tmp_r1 if 'bowtie' in x or 'minimap' in x):
                bowtie = 1
            continue
        # if read is paired to the next file
        if tmp_r1[6] == '=':
            line = next(loop_read)
            tmp_r2 = line.strip().split('\t')
            start += len(line)
            if start > end:
                break
            if tmp_r1[0] == tmp_r2[0]:  # paired reads have the same name
                # identify forward strand
                if int(f'{int(tmp_r1[1]):012b}'[::-1][4]):
                    r1 = tmp_r2
                    r2 = tmp_r1
                else:
                    r1 = tmp_r1
                    r2 = tmp_r2
                read_list = [r1, r2]
            else:
                raise Exception('Paired read not found')
        else:  # only if there is one read
            read_list = [tmp_r1]
        tmpcounts = set()  # initialize name
        read_name = tmp_r1[0]
        indel_size = 0
        tmp_wt_count = [0]*int(len(wt)/3+1)
        read_length = [100000, 0]
        for r3 in read_list:
            if (bin(int(r3[1]))[::-1][1] or bin(int(r3[1]))[::-1][3]) and int(r3[4]) > quality:  # if forward read is aligned and higher than quality threshold
                [tmpcounts, read_length, tmp_wt_count, indel_size] = parse_cigar(r3, wt, quality_nt, tmpcounts, read_length, tmp_wt_count, indel_size, bowtie)
                if len(tmpcounts) > 0:  # only record if mutations found
                    # sort data by position
                    # turn tmpcounts into a list without affecting tmpcounts
                    tmpcounts_too = list(tmpcounts.copy())
                    mutation = [tmpcounts_too[i[0]] for i in sorted(enumerate([int(x.split('_')[0]) for x in tmpcounts_too]), key=lambda x:x[1])]
                    file.write(read_name + ',' + ','.join(mutation) + '\n')
                    # only record if no indels were found
                    if app.indel or (all(['ins' not in x for x in mutation]) and all(['del' not in x for x in mutation])):
                        if app.muts:
                            for mut in mutation:
                                if mut.strip('*') in app.muts_list:
                                    try:
                                        mutation_dict[mut] += 1  # simple dictionary of mutations
                                    except KeyError:
                                        mutation_dict[mut] = 1
                        else:
                            for mut in mutation:
                                try:
                                    mutation_dict[mut.strip('*')] += 1  # simple dictionary of mutations
                                except KeyError:
                                    mutation_dict[mut.strip('*')] = 1
                    if variant_check or correlation_check:
                        if not app.variantfull or (app.variantfull and read_length[0] == 0 and read_length[1] >= len(wt) and indel_size < app.max_indel):
                            # convert to AA
                            variant_mut = []
                            for mut in mutation:
                                if app.muts:
                                    if mut.strip('*') in app.muts_list:
                                        variant_mut.append(mut)
                                else:
                                    variant_mut.append(mut)
                            if variant_mut:
                                # find other mutations not designed
                                other_mut = [x for x in mutation if x not in variant_mut]
                                if '+'.join(variant_mut) in variants.keys():
                                    variants['+'.join(variant_mut)] += 1
                                    # store other_mut in a dictionary
                                    variant_tracking['+'.join(variant_mut)]['read_names'].extend([read_position])
                                    for muto in other_mut:
                                        if muto in variant_tracking['+'.join(variant_mut)].keys():
                                            variant_tracking['+'.join(variant_mut)][muto] += 1
                                        else:
                                            variant_tracking['+'.join(variant_mut)][muto] = 1
                                else:
                                    variants['+'.join(variant_mut)] = 1
                                    # store other_mut in a dictionary
                                    variant_tracking['+'.join(variant_mut)] = {}
                                    variant_tracking['+'.join(variant_mut)]['read_names'] = [read_position]
                                    for muto in other_mut:
                                        variant_tracking['+'.join(variant_mut)][muto] = 1
                            else:
                                for muto in mutation:
                                    if muto in variant_tracking['WT'].keys():
                                        variant_tracking['WT'][muto] += 1
                                        variant_tracking['WT']['read_names'].extend([read_position])
                                    else:
                                        variant_tracking['WT'][muto] = 1
                                        variant_tracking['WT']['read_names'].extend([read_position])
                else:
                    if variant_check and (not app.variantfull or (app.variantfull and read_length[0] == 0 and read_length[1] >= len(wt) and indel_size < app.max_indel)):
                        variants['WT'] += 1
                if sum(tmp_wt_count) > 1:  # and not indel:
                    wt_positions = [xi for xi, x in enumerate(tmp_wt_count) if x]  # all positions wt codon was found
                    for wt_idx in wt_positions:
                        wt_count[wt_idx] += 1  # only add to total count after both reads have been processed
                    if correlation_check:
                        wt_file.write(','.join(str(x) for x in wt_positions) + '\n')
    loop_read.close()
    file.close()
    wt_file.close()
    # os.remove(sam_file)
    # loop through variant_tracking and find mutations that are above a threshold
    if variant_check and app.muts:
        print('Rechecking variants for unknown mutations')
        variants = {'WT': 0}
        app.variant_count = 10
        app.variant_percent = 0.33
        for variant in variant_tracking.keys():
            # total number of reads
            total_reads = len(variant_tracking[variant]['read_names'])
            # find mutations in the dictionary that are above threshold
            if variant == 'WT':
                other_mutations = [k for k, v in variant_tracking[variant].items() if k != 'read_names' and (
                            v > app.variant_count and v / total_reads > 0.01)]
            else:
                other_mutations = [k for k, v in variant_tracking[variant].items() if k != 'read_names' and (v > app.variant_count and v / total_reads > app.variant_percent)]
            # rerun analysis for the reads associated with that mutation
            loop_read = open(sam_file, 'r')
            for read in variant_tracking[variant]['read_names']:
                # go to that read in the file
                loop_read.seek(read)
                line = next(loop_read)
                tmp_r1 = line.split('\t')
                # if read is paired to the next file
                if tmp_r1[6] == '=':
                    line = next(loop_read)
                    tmp_r2 = line.split('\t')
                    start += len(line)
                    if start > end:
                        break
                    if tmp_r1[0] == tmp_r2[0]:  # paired reads have the same name
                        # identify forward strand
                        if int(f'{int(tmp_r1[1]):012b}'[::-1][4]):
                            r1 = tmp_r2
                            r2 = tmp_r1
                        else:
                            r1 = tmp_r1
                            r2 = tmp_r2
                        read_list = [r1, r2]
                    else:
                        raise Exception('Paired read not found')
                else:  # only if there is one read
                    read_list = [tmp_r1]
                tmpcounts = set()  # initialize name
                read_name = tmp_r1[0]
                indel_size = 0
                tmp_wt_count = [0] * int(len(wt) / 3 + 1)
                read_length = [100000, 0]
                for r3 in read_list:
                    [tmpcounts, read_length, tmp_wt_count, indel_size] = parse_cigar(r3, wt, quality_nt,
                                                                                         tmpcounts, read_length,
                                                                                         tmp_wt_count, indel_size,
                                                                                         bowtie)
                    if len(tmpcounts) > 0:  # only record if mutations found
                        # sort data by position
                        tmpcounts = list(tmpcounts)
                        mutation = [tmpcounts[i[0]] for i in
                                    sorted(enumerate([int(x.split('_')[0]) for x in tmpcounts]),
                                           key=lambda x: x[1])]
                        if app.indel or (all(['ins' not in x for x in mutation]) and all(
                                ['del' not in x for x in mutation])):  # only recording if no indels were found
                            if app.muts:
                                for mut in mutation:
                                    if mut.strip('*') in other_mutations:
                                        try:
                                            mutation_dict[mut] += 1  # simple dictionary of mutations
                                        except KeyError:
                                            mutation_dict[mut] = 1
                        if not app.variantfull or (app.variantfull and read_length[0] == 0 and read_length[1] >= len(wt) and indel_size < app.max_indel):
                            # convert to AA
                            variant_mut = []
                            for mut in mutation:
                                if app.muts:
                                    if mut.strip('*') in app.muts_list + other_mutations:
                                        if mut.strip('*').split('_')[1] in codon_chart.keys():
                                            variant_mut.append(mut.split('_')[0] + '_' + codon_chart[mut.strip('*').split('_')[1]])
                                        else:
                                            variant_mut.append(mut)
                            if variant_mut:
                                if '+'.join(variant_mut) in variants.keys():
                                    variants['+'.join(variant_mut)] += 1
                                else:
                                    variants['+'.join(variant_mut)] = 1
                            else:
                                variants['WT'] += 1
    return mutation_dict, variants, wt_count


def write_output(mutation_dict, variants, wt_count, wt, outfile):
    # plot log histogram of mutations
    # plt.figure()
    # plt.hist(mutation_dict.values(), bins=range(1, max(mutation_dict.values()) + 2))
    # plt.xlabel('Number of reads')
    # plt.ylabel('Number of mutations')
    # plt.yscale('log', nonposy='clip')
    # plt.savefig(outfile + '_mutations_hist.png')
    # plt.close()
    # add wt count to mutation dictionary
    for i, wt_c in enumerate(wt_count[1:]):
        mutation_dict[str(i+1)+'_'+str(wt.seq[i*3:i*3+3])] = wt_c
    with open(outfile + '_mutlist.csv', 'w') as f:
        f.write('position,codon,AA,count\n')
        for mut in mutation_dict.keys():
            f.write(','.join(mut.split('_'))+',')
            if mut.split('_')[1] in codon_chart.keys():
                f.write(codon_chart[mut.split('_')[1]]+',')
            elif '*' in mut.split('_')[1]:
                f.write(codon_chart[mut.split('_')[1].strip('*')]+'*,')
            else:
                f.write(mut.split('_')[1]+',')
            f.write(str(mutation_dict[mut])+'\n')
    # remove in and del from mutation_dict
    mutation_dict = {k: v for k, v in mutation_dict.items() if 'ins' not in k and 'del' not in k}
    mutations = pd.DataFrame({'name': list(mutation_dict.keys()), 'count': list(mutation_dict.values())})
    mutations[['position', 'AA']] = mutations['name'].str.split('_', expand=True)
    for row in range(len(mutations)):
        if '*' in mutations.loc[row, 'AA']:
            mutations.loc[row, 'AA'] = codon_chart[mutations.loc[row, 'AA'].strip('*')] + '*'
        else:
            mutations.loc[row, 'AA'] = codon_chart[mutations.loc[row, 'AA']]
    mutations['position'] = pd.to_numeric(mutations['position'])
    matrix = mutations.pivot_table(index='AA', columns='position', values='count', aggfunc='sum')
    matrix.to_csv(outfile + '_matrix.csv')
    if variants:
        # plt.figure()
        # tmp_data = [v for k, v in variants.items() if v > 1 and k != 'WT']
        # plt.hist(tmp_data, bins=range(1, max(tmp_data) + 2))
        # plt.xlabel('Number of reads')
        # plt.ylabel('Number of variants')
        # plt.yscale('log', nonposy='clip')
        # plt.savefig(outfile + '_variants_hist.png')
        # plt.close()
        # Add wt to mutation_dict
        with open(outfile + '_variants.csv', 'w') as f:
            f.write('count,mutation\n')
            for mut in variants.keys():
                f.write(str(variants[mut]) + ',')
                f.write(mut + '\n')


def david_paired_analysis(mut_files, wt_files, app, log):
    #percentage_reads = list(range(0, app.reads, int(app.reads/99)))
    outfile = os.path.splitext(mut_files)[0].split('.fastq')[0]
    ## Mutation Pairs ##
    mutation_pair_dict = {}
    correlations_wt = {}
    single_muts = []
    with open(mut_files, 'r') as mutation_file:
        tmp = next(mutation_file)
        tmp = next(mutation_file)
        for mutation_list in mutation_file:
            mutations = mutation_list.strip('\n').split(',')[1].split('+')
            if len(mutations) > 1:
                for mut_1, mut_2 in itertools.combinations(mutations, 2):
                    # add to total counts
                    mut_ints = sorted([int(mut_1.split('_')[0]), int(mut_2.split('_')[0])])
                    muts = sorted([mut_1, mut_2])
                    try:
                        correlations_wt[','.join(str(intx) for intx in mut_ints)] += int(mutation_list.strip('\n').split(',')[0])
                    except KeyError:
                        correlations_wt[','.join(str(intx) for intx in mut_ints)] = int(mutation_list.strip('\n').split(',')[0])
                    try:
                        mutation_pair_dict['+'.join(muts)] += int(mutation_list.strip('\n').split(',')[0])  # simple dictionary of mutations
                    except KeyError:
                        mutation_pair_dict['+'.join(muts)] = int(mutation_list.strip('\n').split(',')[0])
            #else:
            #    # save single mutations for later
            #    single_muts.append(mutation_list)
    #threshold = 1
    #mutation_pair_dict = {k:v for k,v in mutation_pair_dict.items() if v>threshold}
    #correlations_wt = {k:v for k,v in correlations_wt.items() if v>threshold}

    ## WT pairs ##
    ## used this resource https://nurdabolatov.com/parallel-processing-large-file-in-python
    # Maximum number of processes we can run at a time
    cpu_count = mp.cpu_count()
    file_size = os.path.getsize(wt_files)
    chunk_size = file_size // cpu_count
    chunk_args = chunkit(wt_files, file_size, chunk_size)
    # parallel
    with mp.Pool(cpu_count) as pool:
        answers = pool.starmap(wt_count, chunk_args)
    for answer in answers:
        correlations_wt = {k: correlations_wt.get(k, 0) + answer.get(k, 0) for k in set(correlations_wt) | set(answer)}
    wt_matrix = pd.DataFrame({'pair': list(correlations_wt.keys()), 'count': list(correlations_wt.values())})
    wt_matrix[['pos1', 'pos2']] = wt_matrix['pair'].str.split(',', expand=True)
    wt_matrix['pos1'] = pd.to_numeric(wt_matrix['pos1'])
    wt_matrix['pos2'] = pd.to_numeric(wt_matrix['pos2'])
    matrix = wt_matrix.pivot_table(index='pos1', columns='pos2', values='count', aggfunc=sum)
    with open(outfile+'_mut_pairs.csv', 'w') as f:
        f.write('mutants,count,percentage\n')
        for mut in mutation_pair_dict.keys():
            f.write(mut+',')
            f.write(str(mutation_pair_dict[mut])+',')
            # calculate percentage of reads for each pair
            mut_1, mut_2 = mut.split('+')
            try:
                f.write(str((mutation_pair_dict[mut]+1) / (correlations_wt[mut_1.split('_')[0]+','+mut_2.split('_')[0]]+1))+'\n')
            except KeyError:
                f.write('\n')


def wt_count(file_name, chunk_start, chunk_end):
    correlations_wt = {}
    with open(file_name, 'r') as f:
        f.seek(chunk_start)
        for positions in f:
            chunk_start += len(positions)
            if chunk_start > chunk_end:
                break
            for wt_pos in itertools.combinations(positions.strip('\n').split(',')[1:], 2):
                try:
                    correlations_wt[','.join(wt_pos)] += 1
                except KeyError:
                    correlations_wt[','.join(wt_pos)] = 1
    return correlations_wt



def align_all_bbmap(sequencing_file, reference, sam_file, par, bbmap_script=f'java -cp {script_path} align2.BBMap',
                    max_gap=500, paired_sequencing_file=False):
    if paired_sequencing_file:
        if par:
            command = f"{bbmap_script} ref={reference} in={sequencing_file} " \
                      f"in2={paired_sequencing_file} maxindel={max_gap} t={par} local outm={sam_file}"
        else:
            command = f"{bbmap_script} ref={reference} in={sequencing_file} " \
                  f"in2={paired_sequencing_file} maxindel={max_gap} local outm={sam_file}"
    else:
        if par:
            command = f"{bbmap_script} ref={reference} in={sequencing_file} " \
                      f"maxindel={max_gap} local t={par} outm={sam_file}"
        else:
            command = f"{bbmap_script} ref={reference} in={sequencing_file} " \
                    f"maxindel={max_gap} local outm={sam_file}"
    print(command)
    ret = os.system(command)
    assert ret == 0

def align_long_read(sequencing_file, reference, sam_file, par, script=minimap_path):
    command = f"{script} {reference} {sequencing_file} -x map-ont " \
              f"-O 15,24 -a --eqx > {sam_file}"
    print(command)
    ret = os.system(command)
    assert ret == 0


def align_domain(app, sequencing_file, reference, domains_file, domain_out_file):
    # read domains
    domains = list(Bio.SeqIO.parse(domains_file, 'fasta'))
    # generate fragments small enough to find ends but large enough to distinguish from other domains
    for i, domain in enumerate(domains):
        end5 = 19
        prime5 = domain.seq[:end5]
        while any([prime5 in x.seq for x in domains[:i]+domains[i+1:]]):
            end5 += 1
            prime5 = domain.seq[:end5]
        end3 = -20
        prime3 = domain.seq[end3:]
        while any([prime3 in x.seq for x in domains[:i]+domains[i+1:]]):
            end3 -= 1
            prime3 = domain.seq[end3:]
    domain_counts = pd.DataFrame()
    quality = int(app.minq)
    quality_nt = int(app.minb)
    wt = Bio.SeqIO.read(reference, 'fasta').seq
    wt_count = [0] * int(len(wt) / 3 + 1)
    fraction = 100/len(domains)
    for domain in domains:
        domain_reference = os.path.join(os.path.dirname(domains_file), domain.id + '.fasta')
        # create a reference file for each domain
        with open(domain_reference, 'w') as outfile:
            outfile.write(f">{domain.id}\n{domain.seq}\n")
            outfile.close()
        print('##### Retaining reads that match insert #####')
        insert_filtered_fnames = os.path.join(os.path.dirname(domains_file), domain.id + '_insert_filtered.fastq')
        command = [f'java -cp {script_path} jgi.BBDuk',
                   'ref=%s' % domain_reference,
                   'in=%s' % sequencing_file,
                   'outm=%s' % insert_filtered_fnames,
                   'mink=%d' % end5,
                   'mm=f kmask=N']
        print('Running: \n\t%s' % ' '.join(command))
        ret = os.system(' '.join(command))
        assert ret == 0, 'bbduk.sh failed for %s!' % domain_reference

        print('##### Retaining reads also matching backbone #####')
        insert_bbone_filtered_fnames = os.path.join(os.path.dirname(domains_file), domain.id + '_insert_bbone_filtered.fastq')
        command = [f'java -cp {script_path} jgi.BBDuk',
                   'ref=%s' % reference,
                   'in=%s' % insert_filtered_fnames,
                   'outm=%s' % insert_bbone_filtered_fnames,
                   'mink=%d' % 12,
                   'hdist=1']
        print('Running: \n\t%s' % ' '.join(command))
        ret = os.system(' '.join(command))
        assert ret == 0, 'bbduk.sh failed for %s!' % reference

        print('##### Find Domain position in read #####')
        # if long read use minimap2

        filtered_masked_fnames = os.path.join(os.path.dirname(domains_file), domain.id + '_domain.sam')
        command = [f'java -cp {script_path} align2.BBMap',
                   'ref=%s' % domain_reference,
                   'in=%s' % insert_bbone_filtered_fnames,
                   'out=%s' % filtered_masked_fnames,
                   'touppercase=t nodisk overwrite local vslow minratio=0.01 k=8']
        ret = os.system(' '.join(command))
        assert ret == 0, 'bbmap.sh failed for %s!' % domain_reference

        # write to csv and create table
        read_length = [100000, 0]
        file_out = open(domain_reference.split('.fasta')[0] + '.csv', 'w')
        with open(filtered_masked_fnames, 'r') as file:
            for row in file:
                if '@' in row[0]:
                    continue
                # if read is paired to the next file
                tmp_r1 = row.split('\t')
                if tmp_r1[6] == '=':
                    tmp_r2 = next(file).split('\t')
                    if tmp_r1[0] == tmp_r2[0]:  # paired reads have the same name
                        read_list = [tmp_r1, tmp_r2]
                    else:
                        raise Exception('Paired read not found')
                else:  # only if there is one read
                    read_list = [tmp_r1]
                tmpcounts = {tmp_r1[0]}  # initialize name
                for r3 in read_list:
                    if (bin(int(r3[1]))[::-1][1] or bin(int(r3[1]))[::-1][3]) and '*' not in r3[5]: # make sure it aligned
                        if all([tmp == 'N' for tmp in r3[9][r3[9].find('N'):r3[9].rfind('N')]]):  # make sure no mutations or gaps
                            # set index
                            ref = int(r3[3]) - 1  # subtract 1 for python numbering
                            read_length[0] = min(read_length[0], ref)
                            read = 0
                            # cycle through each codon until end of read or end of reference sequence
                            fragment_cigar = list(filter(None, re.split(r'(\d+)', r3[5])))
                            codon_start = ref + ref % 3
                            codon_end = ref
                            for cigar in [fragment_cigar[ii:ii + 2] for ii in range(0, len(fragment_cigar), 2)]:
                                codon_end += int(cigar[0])
                                if 'X' == cigar[1]:
                                    # this is a substitution
                                    ref += int(cigar[0])
                                    read += int(cigar[0])
                                elif '=' == cigar[1]:
                                    # this is a perfect match
                                    ref += int(cigar[0])
                                    read += int(cigar[0])
                                elif 'D' == cigar[1]:
                                    # this is a deletion
                                    # insert '-' for gap
                                    r3[9] = r3[9][:read] + '-' * int(cigar[0]) + r3[9][read:]
                                    r3[10] = r3[10][:read] + '-' * int(cigar[0]) + r3[9][read:]
                                    read += int(cigar[0])
                                    ref += int(cigar[0])
                                elif 'I' == cigar[1]:
                                    # this is an insertion
                                    read += int(cigar[0])
                                elif 'S' == cigar[1]:
                                    # soft clipping
                                    read += int(cigar[0])
                                elif 'M' == cigar[1]:
                                    # bad read (N) - this is the domain
                                    if int(cigar[0]) > end5:
                                        domain_position = int(ref / 3) + 1
                                        if read == 0 or int(cigar[0]) > len(domain.seq):
                                            domain_position = int((ref + int(cigar[0]) - len(domain.seq))/3) + 1
                                        try:
                                            domain_counts.loc[domain_position, domain.id] += 1
                                        except KeyError:
                                            domain_counts.loc[domain_position, domain.id] = 1
                                        tmpcounts.add(str(domain_position))
                                    ref += int(cigar[0])
                                    read += int(cigar[0])
                                else:
                                    raise TypeError('Cigar format not found')
                            if codon_end - codon_start >= 3:
                                # make sure this doesn't exceed the length of the gene using min
                                wt_range = range(int(codon_start / 3) + 1, min(int(codon_end/3) + 1, int(len(wt) / 3)), 1)  # need to add 1 for range
                                for wt_i in wt_range:
                                    wt_count[wt_i] += 1
                if len(tmpcounts) > 1:  # only record if mutations found
                    # append data
                    sorted_counts = sorted(tmpcounts, reverse=True)
                    file_out.write(','.join(sorted_counts) + '\n')
        #os.remove(aligned_fnames)
        #os.remove(insert_filtered_fnames)
        #os.remove(insert_bbone_filtered_fnames)
        #os.remove(filtered_masked_fnames)
        #os.remove(domain_reference)
    for wt_pos, count in enumerate(wt_count):
        domain_counts.loc[wt_pos, 'wt'] = count
    # convert NaN to 0
    domain_counts = domain_counts.fillna(0)
    domain_counts.to_csv(domain_out_file.split('.sam')[0]+'_table.csv')
    file_out.close()
    return domain_counts


def map_structure(structure_file, matrix, min_q, max_q, average, name):
    # pymol.finish_launching()
    pymol.cmd.load(structure_file)
    structure_name = os.path.basename(structure_file).split('.')[0]
    pymol.stored.list = []
    pymol.cmd.iterate("(name ca)", "stored.list.append(resi)")
    matrix2 = matrix.iloc[:, 2:].reindex(columns=[int(x) for x in pymol.stored.list], fill_value=0).fillna(0)
    #print([float(x) for x in list(matrix2.iloc[22, :])])
    pymol.stored.all = [float(x) for x in list(matrix2.iloc[22, :])]  # average row
    pymol.cmd.alter(selection=structure_name, expression='b=0.0')
    pymol.cmd.alter(selection=structure_name + ' and n. CA', expression='b=stored.all.pop(0)')
    pymol.cmd.spectrum(expression='b', palette='red_white_green', selection=structure_name + ' and n. CA',
                       minimum=min_q, maximum=max_q)
    pymol.cmd.save(name + '_mapping.pse')
    # loop through amino acids groups and store average in pymol b factor
    aa_groups = {'Polar': ['N', 'Q', 'S', 'T'], 'Hydrophobic': ['A', 'I', 'L', 'M', 'V'], 'Positive': ['R', 'K', 'H'], 'Negative': ['D', 'E'], 'Aromatic': ['F', 'Y', 'W']}
    for aa_group, aa_list in aa_groups.items():
        sub_matrix = matrix2.loc[aa_list, :]
        # replace 0 in sub_matrix with numpy NaN
        sub_matrix = sub_matrix.replace(0, np.nan)
        # calculate average and ignore NaN
        aa_average = [float(x) for x in list(sub_matrix.mean(axis=0, skipna=True))]
        # replace NaN with average
        aa_average = [average if math.isnan(x) else x for x in aa_average]  # this will maintain 0 as white
        pymol.stored.all = aa_average
        pymol.cmd.alter(selection=structure_name, expression='b=0.0')
        pymol.cmd.alter(selection=structure_name + ' and n. CA', expression='b=stored.all.pop(0)')
        pymol.cmd.spectrum(expression='b', palette='red_white_green', selection=structure_name + ' and n. CA',
                           minimum=min_q, maximum=max_q)
        pymol.cmd.save(name + '_mapping_' + aa_group + '.pse')
        # select enriched residues for given aa group
        enriched = sub_matrix.loc[:, sub_matrix.mean(axis=0, skipna=True) > average]
        pymol.cmd.select('enriched_' + aa_group, structure_name + ' and n. CA and resi ' + '+'.join([str(x) for x in list(enriched.columns)]))
        pymol.cmd.set('dot_solvent', 1)
        surface_area = pymol.cmd.get_area('enriched_' + aa_group, state=1, load_b=0)/len(enriched.columns)
        # print surface area
        print('Surface area of ' + aa_group + ' enriched residues: ' + str(surface_area))
        # select enriched residues over 80% percentile for given aa group
        enriched = sub_matrix.loc[:, sub_matrix.mean(axis=0, skipna=True) > np.nanpercentile(sub_matrix, 80)]
        pymol.cmd.select('enriched_80_' + aa_group, structure_name + ' and n. CA and resi ' + '+'.join([str(x) for x in list(enriched.columns)]))
        surface_area = pymol.cmd.get_area(selection='enriched_80_' + aa_group, state=1, load_b=0)/len(enriched.columns)
        # print surface area
        print('Surface area of ' + aa_group + ' enriched residues over 80% percentile: ' + str(surface_area))



def calc_enrichment(app, root, base, selected, name, mincount, wt_file, structure_file, muts_file):
    wt_nt_seq = Bio.SeqIO.read(wt_file, "fasta").upper()
    # translate wt_nt_seq to AA into a dataframe
    wt_seq = pd.DataFrame({'WT_AA': wt_nt_seq.translate()}).T
    if muts_file:
        # filter base for only codons in muts_file
        muts = pd.read_csv(muts_file, sep='\t', header=None)
        # split text into position and codon
        muts = muts[0].str.split('_', expand=True)
        muts.columns = ['position', 'codon']
        # add wt codons to muts
        for idx, i in enumerate(range(1, len(wt_nt_seq), 3)):
            muts = pd.concat([muts, pd.DataFrame({"position":[idx+1], "codon":[str(wt_nt_seq.seq[i-1:i+2])]})], ignore_index=True)
        muts['position'] = muts['position'].astype(int)
        muts['codon'] = muts['codon'].astype(str)
        # filter base and selected for only positions and codons in muts
        # must match the same position and codon in the same row
        base = base.merge(muts, on=['position', 'codon'])
        selected = selected.merge(muts, on=['position', 'codon'])
        base = base.groupby([x for x in base.columns if x not in ['count', 'codon']], as_index=False).sum()
        selected = selected.groupby([x for x in selected.columns if x not in ['count', 'codon']], as_index=False).sum()
    else:
        if not selected.empty:
            base = base.groupby([x for x in base.columns if x not in ['count', 'codon']], as_index=False).agg('sum')
            selected = selected.groupby([x for x in selected.columns if x not in ['count', 'codon']], as_index=False)['count'].agg('sum')
        else:
            # loop through base and groupby position
            tmp_base = base[0]
            tmp_base = tmp_base.groupby([x for x in tmp_base.columns if x not in ['count', 'codon']], as_index=False)['count'].agg('sum')
            df = tmp_base[[x for x in tmp_base.columns if x not in ['codon']]]
            for tmp_base in base[1:]:
                tmp_base = tmp_base.groupby([x for x in tmp_base.columns if x not in ['count', 'codon']], as_index=False)['count'].agg('sum')
                df = pd.merge(df, tmp_base[[x for x in tmp_base.columns if x not in ['codon']]], left_on=['position', 'AA'], right_on=['position', 'AA'], how='outer')
    read_count = 0
    # combine into one dataframe
    if 'position' in list(base.columns):
        if not selected.empty:
            df = pd.merge(base, selected, left_on=['position', 'AA'], right_on=['position', 'AA'], suffixes=('_base', '_select'), how='outer')
            df = df.fillna(0)
            percentage_reads = list(range(0, len(df), int(len(df) / 99) + (len(df) % 99 > 0)))
            root.update_idletasks()
            app.progress['value'] = 0
            root.update()
            # by AA position
            for i in df.position.unique():
                selection = df.loc[df.position == i, :]
                inp_count = sum(selection.loc[:, 'count_base'])+1
                sel_count = sum(selection.loc[:, 'count_select'])+1
                for select, row in selection.iterrows():
                    inp_score = selection.loc[select, 'count_base'] + 1
                    sel_score = selection.loc[select, 'count_select'] + 1
                    score = math.log(sel_score / sel_count) - math.log(inp_score / inp_count)
                    se = math.sqrt((1 / sel_score) + (1 / inp_score) + (1 / sel_count) + (1 / inp_count))
                    df.loc[select, 'base_frequency'] = inp_score / inp_count
                    df.loc[select, 'select_frequency'] = sel_score / sel_count
                    df.loc[select, 'score'] = score
                    df.loc[select, 'std_error'] = se
                    read_count += 1
                    if read_count in percentage_reads:
                        root.update_idletasks()
                        app.progress['value'] += 1
                        root.update()
        else:
            tmp_df = pd.DataFrame()
            linear_regressor = LinearRegression()
            for i in df.position.unique():
                selection = df.loc[df.position == i, :]
                # loop through each count and calculate frequency
                for column in [x for x in selection.columns if x not in ['position', 'AA']]:
                    selection.loc[:, column+'_freq'] = selection[column] / sum(selection[column])
                # calculate slope of each row with the freq columns
                freqs = selection[[x for x in selection.columns if '_freq' in x]]
                for irow in freqs.index:
                    x = np.arange(0, len(freqs.columns))
                    y = freqs.loc[irow, ].values
                    # find index of nan in y
                    nan_index = np.argwhere(np.isnan(y))
                    if len(nan_index) < len(y) - 1:
                        # remove nan from x and y
                        x = np.delete(x, nan_index)
                        y = np.delete(y, nan_index)
                        fit = linear_regressor.fit(x.reshape(-1, 1), y.reshape(-1, 1))
                        tmp_df.loc[irow, 'score'] = fit.coef_[0][0]
            df = pd.merge(df, tmp_df, left_index=True, right_index=True)
        df.to_csv(name, index=False)
        df['mutation'] = df['position'].astype(str) + '_' + df['AA']
        # create matrix
        if not selected.empty:
            df = df.loc[(df['count_base'] > mincount) | (df['count_select'] > mincount), ]
        else:
            df = df.loc[df[[x for x in df.columns if 'count' in x]].sum(axis=1) > mincount, ]
        matrix = df.pivot_table(index='AA', columns='position', values='score', aggfunc=np.mean)

        # color matrix based on count_matrix
        matrix = pd.concat([wt_seq, matrix])

        for i in matrix.keys().values.tolist():
            matrix.loc[matrix.loc['WT_AA', i], i] = np.nan
        # synonymous mutations replacement
        for row in [x for x in matrix.index if not pd.isnull(x) and '*' in x and not '*' == x]:
            for col in matrix.columns:
                if not np.isnan(matrix.loc[row, col]):
                    matrix.loc[row.strip('*'), col] = matrix.loc[row, col]
        aa_order = ['WT_AA', 'A', 'V', 'I', 'L', 'M', 'F', 'Y', 'W', 'R', 'H', 'K', 'D', 'E', 'S', 'T', 'N', 'Q', 'C', 'G', 'P', '*']
        matrix = matrix.reindex(aa_order)

        matrix['Average'] = matrix.drop('WT_AA').mean(axis=1)
        # move average to the beginning
        matrix = matrix[['Average'] + [c for c in matrix if c not in ['Average']]]
        matrix.loc['AveragePosition'] = matrix.drop(['*', 'WT_AA']).mean(axis=0)
        # loop through amino acid groups and calculate average
        aa_groups = {'Polar': ['N', 'Q', 'S', 'T'], 'Hydrophobic': ['A', 'I', 'L', 'M', 'V'], 'Aromatic': ['F', 'Y', 'W'],
                     'Positive': ['R', 'K', 'H'], 'Negative': ['D', 'E']}
        for aa_group in aa_groups:
            matrix.loc[aa_group] = matrix.loc[aa_groups[aa_group]].mean(axis=0)
        writer = pd.ExcelWriter(name.replace('results.csv', '') + 'combined_matrix.xlsx', engine='xlsxwriter')
        matrix.to_excel(writer, sheet_name='Enrichment')
        if not selected.empty:
            count_matrix = df.pivot_table(index='AA', columns='position',
                                          values='count_base', aggfunc=np.sum) + \
                           df.pivot_table(index='AA', columns='position',
                                          values='count_select', aggfunc=np.sum)
            count_matrix = pd.concat([count_matrix,wt_seq])
            count_matrix = count_matrix.reindex(aa_order)
            count_matrix.to_excel(writer, sheet_name='Counts')
        writer.close()
        # align to structure
        if structure_file:
            quartile = matrix.loc['AveragePosition'].quantile([0.10, 0.90])
            average = matrix.loc['AveragePosition'].mean()
            map_structure(structure_file, matrix, quartile[0.1], quartile[0.9], average, name)
    if 'position' not in list(base.columns):
        # by variant
        df = pd.merge(base, selected, left_on=['mutation'], right_on=['mutation'],
                      suffixes=('_base', '_select'), how='outer')
        df = df.fillna(0)
        percentage_reads = list(range(0, len(df), int(len(df) / 99) + (len(df) % 99 > 0)))
        root.update_idletasks()
        app.progress['value'] = 0
        root.update()
        inp_count = sum(df.loc[:, 'count_base'])+1
        sel_count = sum(df.loc[:, 'count_select'])+1
        for select, row in df.iterrows():
            inp_score = df.loc[select, 'count_base'] + 1
            sel_score = df.loc[select, 'count_select'] + 1
            score = math.log(sel_score / sel_count) - math.log(inp_score / inp_count)
            se = math.sqrt((1 / sel_score) + (1 / inp_score) + (1 / sel_count) + (1 / inp_count))
            df.loc[select, 'score'] = score
            df.loc[select, 'std_error'] = se
            read_count += 1
            if read_count in percentage_reads:
                root.update_idletasks()
                app.progress['value'] += 1
                root.update()
        df.to_csv(name, index=False)
    # volcano plot
    p = ggplot(df.loc[(df.loc[:, 'count_base']+df.loc[:, 'count_select']) > mincount, :]) + aes(x='score', y='std_error') + geom_point() + theme_minimal() + geom_text(
        data=df.loc[df['score'] > 2]) + aes(x='score', y='std_error', label='mutation') + scale_y_reverse()
    print(p)
    print('Finished enrichment calculations')


def combine(app, root, dataset, combineBy, name, threshold, structure_file):
    df = dataset[0]
    if 'position' in list(dataset[0].columns):
        for data in dataset[1:]:
            df = pd.merge(df, data, on=['position', 'AA'], how='outer')
        df['mutation'] = df['position'].astype(str) + '_' + df['AA']
    else:
        for data in dataset[1:]:
            df = pd.merge(df, data, on='mutation', how='outer')
    # correlation plot first
    if combineBy:
        tmpdf = df.loc[(df[[x for x in df.columns if 'std_error' in x]] < threshold).apply(sum, axis=1) > len([x for x in df.columns if 'std_error' in x])/2]
    else:
        tmpdf = df.loc[(df[[x for x in df.columns if 'count' in x]] > threshold).apply(sum, axis=1) > len([x for x in df.columns if 'count' in x])/2]
    g = sns.PairGrid(tmpdf[[x for x in df.columns if 'score' in x]])
    g.map_diag(sns.histplot)
    g.map_offdiag(sns.regplot)
    # add r2 value
    for i, j in zip(*np.triu_indices_from(g.axes, 1)):
        g.axes[i, j].annotate('r2 = ' + str(round(tmpdf[[x for x in df.columns if 'score' in x]].corr().iloc[i, j] ** 2, 3)), (0.1, 0.9), xycoords='axes fraction', ha='center', va='center')
    plt.show()
    plt.pause(1)
    read_count = 0
    percentage_reads = list(range(0, len(df), int(len(df) / 99) + (len(df) % 99 > 0)))
    root.update_idletasks()
    app.progress['value'] = 0
    root.update()
    # sum count_base
    df['count_combined_base'] = df[[x for x in df.columns if 'count_base' in x]].sum(axis=1, min_count=1)
    # sum count_select
    df['count_combined_select'] = df[[x for x in df.columns if 'count_select' in x]].sum(axis=1, min_count=1)
    if not combineBy:
        df = df.fillna(0)
        if 'position' in list(df.columns):
            # by AA position
            wt_seq = pd.DataFrame()
            for i in df.position.unique():
                selection = df.loc[df.position == i, :]
                wt_seq.loc['WT_AA', i] = (selection.loc[selection['count_combined_base'].idxmax(), 'AA'])
                inp_count = np.nansum(selection.loc[:, 'count_combined_base'] + 1)
                sel_count = np.nansum(selection.loc[:, 'count_combined_select'] + 1)
                for select, row in selection.iterrows():
                    inp_score = selection.loc[select, 'count_combined_base'] + 1
                    sel_score = selection.loc[select, 'count_combined_select'] + 1
                    score = math.log(sel_score / sel_count) - math.log(inp_score / inp_count)
                    se = math.sqrt((1 / sel_score) + (1 / inp_score) + (1 / sel_count) + (1 / inp_count))
                    if app.multiply.get():
                        df.loc[select, 'combined_score'] = score * (-(1 / (inp_score + sel_score)) ** 0.5 + 1)
                    else:
                        df.loc[select, 'combined_score'] = score
                    df.loc[select, 'combined_std_error'] = se
                    read_count += 1
                    if read_count in percentage_reads:
                        root.update_idletasks()
                        app.progress['value'] += 1
                        root.update()
                df.loc[df.position == i, 'WT_AA'] = selection.loc[selection[[x for x in selection.keys() if 'count_base' in x][0]].idxmax(), 'AA']
        else:
            # by variant
            inp_count = np.nansum(df.loc[:, 'count_combined_base'] + 1)
            sel_count = np.nansum(df.loc[:, 'count_combined_select'] + 1)
            for select, row in df.iterrows():
                inp_score = df.loc[select, 'count_combined_base'] + 1
                sel_score = df.loc[select, 'count_combined_select'] + 1
                score = math.log(sel_score / sel_count) - math.log(inp_score / inp_count)
                se = math.sqrt((1 / sel_score) + (1 / inp_score) + (1 / sel_count) + (1 / inp_count))
                if app.multiply.get():
                    df.loc[select, 'combined_score'] = score * (-(1 / (inp_score + sel_score)) ** 0.5 + 1)
                else:
                    df.loc[select, 'combined_score'] = score
                df.loc[select, 'combined_std_error'] = se
                read_count += 1
                if read_count in percentage_reads:
                    root.update_idletasks()
                    app.progress['value'] += 1
                    root.update()
        df = df.loc[(df['count_combined_base'] > threshold) | (df['count_combined_select'] > threshold), ]
        df.loc[:, 'log_count'] = np.log10(df['count_combined_base'] + df['count_combined_select'])
        # volcano plot
        p = ggplot(df) + aes(x='combined_score', y='log_count') + geom_point() + theme_minimal() + \
            geom_text(data=df.loc[df['combined_score'] > 2]) + aes(x='combined_score', y='log_count', label='mutation')
        print(p)
    else:
        for index, row in df.iterrows():
            scores = row[[x for x in row.index if 'score' in x]]
            se = row[[x for x in row.index if 'std_error' in x]]
            if sum(scores.isna()) < 2 | sum(se.isna()) < 2:
                df.loc[index, 'combined_score'] = scores.mean(skipna=True)
                df.loc[index, 'combined_std_error'] = se.mean(skipna=True)
            else:
                score = list(scores[-scores.isna()])
                se = list(se[-se.isna()])
                answer = rml_estimator(score, se, 500)  # use max 500 estimates to converge
                var_betaML = answer[0]  # return variance maximum likelihood
                eps = answer[1]  # high epsilon means the rml_estimator failed to converge (should filter this)
                if (eps < 0.001):
                    inp_score = df.loc[index, 'count_combined_base'] + 1
                    sel_score = df.loc[index, 'count_combined_select'] + 1
                    score = sum(score * (var_betaML+[x**2 for x in se]) ** -1) / sum((var_betaML+[x**2 for x in se]) ** -1)
                    if app.multiply.get():
                        df.loc[index, 'combined_score'] = score * (-(1 / (inp_score + sel_score)) ** 0.5 + 1)
                    else:
                        df.loc[index, 'combined_score'] = score
                    df.loc[index, 'combined_std_error'] = np.sqrt(1 / sum([x**-2 for x in se]))
                    df.loc[index, 'combined_eps'] = eps
                else:
                    df.loc[index, 'combined_score'] = np.nan
                    df.loc[index, 'combined_std_error'] = np.nan
                    df.loc[index, 'combined_eps'] = np.nan
            read_count += 1
            if read_count in percentage_reads:
                root.update_idletasks()
                app.progress['value'] += 1
                root.update()
        if 'position' in list(df.columns):
            # by AA position
            wt_seq = pd.DataFrame()
            for i in df.position.unique():
                selection = df.loc[df.position == i, :]
                if not np.isnan(selection[[x for x in selection.keys() if 'count_base' in x][0]].idxmax()):
                    wt_seq.loc['WT_AA', i] = (selection.loc[selection[[x for x in selection.keys() if 'count_base' in x][0]].idxmax(), 'AA'])
                    df.loc[df.position == i, 'WT_AA'] = selection.loc[selection[[x for x in selection.keys() if 'count_base' in x][0]].idxmax(), 'AA']
                read_count += 1
                if read_count in percentage_reads:
                    root.update_idletasks()
                    app.progress['value'] += 1
                    root.update()
        df = df.loc[df['combined_std_error'] < threshold, ]
        # volcano plot
        p = ggplot(df) + aes(x='combined_score', y='combined_std_error') + geom_point() + theme_minimal() + \
            geom_text(data=df.loc[df['combined_score'] > 2]) + aes(x='combined_score', y='combined_std_error',
                                                                   label='mutation') + scale_y_reverse()
        print(p)
    df.to_csv(name+'/combined_results.csv', index=False)
    # create matrix
    if 'position' in df.columns:
        matrix = df.pivot_table(index='AA', columns='position', values='combined_score', aggfunc=np.mean)
        # addition of count_combined_base and count_combined_select into matrix
        count_matrix = df.pivot_table(index='AA', columns='position',
                                      values='count_combined_base', aggfunc=np.sum) + \
                       df.pivot_table(index='AA', columns='position',
                                      values='count_combined_select', aggfunc=np.sum)
        missing_positions = list(set(range(1, max(df['position'])+1)).difference(matrix.columns))
        matrix[missing_positions] = np.nan
        matrix = pd.concat([matrix,wt_seq])
        for i in range(1, len(matrix.columns)):
            matrix.loc[matrix.loc['WT_AA', i], i] = np.nan
        # remove nan row
        matrix = matrix.dropna(axis=0, how='all')
        # synonymous mutations replacement
        for row in [x for x in matrix.index]:
            if '*' in row and not '*' == row:
                for col in matrix.columns:
                    if not np.isnan(matrix.loc[row, col]):
                        matrix.loc[row.strip('*'), col] = matrix.loc[row, col]
        aa_order = ['WT_AA', 'A', 'V', 'I', 'L', 'M', 'F', 'Y', 'W', 'R', 'H', 'K', 'D', 'E', 'S', 'T', 'N', 'Q', 'C', 'G', 'P', '*']
        matrix = matrix.reindex(aa_order)
        count_matrix = count_matrix.reindex(aa_order)
        # calculate average
        matrix['Average'] = matrix.drop('WT_AA').mean(axis=1, skipna=True)
        # move average to the beginning
        matrix = matrix[['Average'] + [c for c in matrix if c not in ['Average']]]
        matrix.loc['AveragePosition'] = matrix.drop(['*', 'WT_AA']).mean(axis=0, skipna=True)
        # loop through amino acid groups and calculate average
        aa_groups = {'Polar': ['N', 'Q', 'S', 'T'], 'Hydrophobic': ['A', 'I', 'L', 'M', 'V'], 'Aromatic': ['F', 'Y', 'W'],
                     'Positive': ['R', 'K', 'H'], 'Negative': ['D', 'E']}
        for aa_group in aa_groups:
            matrix.loc[aa_group] = matrix.loc[aa_groups[aa_group]].mean(axis=0, skipna=True)
        writer = pd.ExcelWriter(name + '/combined_matrix.xlsx', engine='xlsxwriter')
        matrix.to_excel(writer, sheet_name='Enrichment')
        count_matrix.to_excel(writer, sheet_name='Counts')
        writer.close()
        # cluster data
        # calculate UMAP and plot
        #umap_data = matrix.drop(['*', 'WT_AA'], axis=0).dropna(axis=1, how='all')
        # remove columns with more than 3 NaNs
        #umap_data = umap_data.loc[:, umap_data.isna().sum(axis=0) < 3]
        # replace nan with 0
        #umap_data = umap_data.fillna(0).transpose()
        #mapper = umap.UMAP().fit(umap_data)
        # plot mapper
        #umap.plot.points(mapper)
        # align to structure
        if structure_file:
            quartile = matrix.loc['AveragePosition'].quantile([0.10, 0.90])
            average = matrix.loc['AveragePosition'].mean()
            map_structure(structure_file, matrix, quartile[0.1], quartile[0.9], average, name)
    print('Finished Combining Scores')



def rml_estimator(y, sigma2i, iterations=500):  # function for recursive maximum likelihood estimator
    iteration = 1
    eps = 100
    sigma2ML = sum((y - np.mean(y)) ** 2) / (len(y) - 1)
    while iteration < iterations:
        if eps < 0.0001 and iteration > 50:  # must get to 50 iterations
            break
        sigma2ML_new = sigma2ML*sum(((y - np.mean(y)) ** 2) * (sigma2ML+sigma2i)**-2) / (sum((sigma2ML+sigma2i)**-1) - (sum((sigma2ML+sigma2i)**-2)/sum((sigma2ML+sigma2i)**-1)))
        eps = abs(sigma2ML - sigma2ML_new)
        sigma2ML = sigma2ML_new
        iteration = iteration+1
        var_betaML = 1 / sum(1 / (sigma2i + sigma2ML))
    return [var_betaML, eps]


def find_barcodes(sequencing_file, paired_sequencing_file, adapters):
    print('##### Determining Barcode Location #####')
    if paired_sequencing_file:
        out1 = sequencing_file.split('.fastq.gz')[0]+'_corrected.fastq.gz'
        command = [f'java -cp {script_path} jgi.BBMerge',
                   'in=%s' % sequencing_file,
                   'in2=%s' % paired_sequencing_file,
                   'out1=%s' % out1,
                   "adapters=%s" % adapters]
        ret = os.system(' '.join(command))
        assert ret == 0, 'finding barcodes failed!'
    # filter reads that align to adapters
    out2 = os.path.join(os.path.dirname(sequencing_file), 'filtered_reads.fastq')
    command = [f'java -cp {script_path} jgi.BBDuk',
               'in=%s' % out1,
               'ref=%s' % adapters,
               'outm=%s' % out2]
    ret = os.system(' '.join(command))
    assert ret == 0, 'Right trimming failed!'
    with open(adapters, 'r') as file:
        next(file)
        left_adapter = next(file).strip()
        next(file)
        right_adapter = next(file).strip()
    tmpbarcodes = os.path.join(os.path.dirname(sequencing_file), 'tmp_barcodes.fastq')
    command = [f'java -cp {script_path} jgi.BBDuk',
               'in=%s' % out2,
               'literal=%s' % left_adapter,
               'out=%s' % tmpbarcodes,
               'ktrim=l k=23 mink=11 hdist=1 tpe tbo']
    ret = os.system(' '.join(command))
    assert ret == 0, 'Right trimming failed!'
    barcodes = os.path.join(sequencing_file.split('.fastq')[0], 'barcodes.fastq')
    command = [f'java -cp {script_path} jgi.BBDuk',
               'in=%s' % tmpbarcodes,
               'literal=%s' % right_adapter,
               'out=%s' % barcodes,
               'ktrim=r k=23 mink=11 hdist=1 tpe tbo trimq=20']
    ret = os.system(' '.join(command))
    assert ret == 0, 'Left trimming failed!'
    ## read barcodes and count
    n = 4
    with open(barcodes, 'r') as fh:
        next(fh)
        line = next(fh)
        count_barcodes = {line.strip(): 1}
        count = 0
        for line in fh:
            count += 1
            if count == n:
                try:
                    count_barcodes[line.strip()] += 1
                except KeyError:
                    count_barcodes[line.strip()] = 1
                count = 0
    with open(os.path.join(sequencing_file.split('.fastq')[0], 'barcode_counts.csv'), 'w') as file:
        for barcode in count_barcodes.keys():
            if count_barcodes[barcode] > 1:
                file.write(barcode+', '+str(count_barcodes[barcode])+'\n')
    print('Finished finding barcodes')


def analyze_sanger(seq_file):
    handle = open(seq_file, 'rb')
    record = SeqIO.read(handle, "abi")
    index = list(record.annotations['abif_raw']['PLOC1'])
    base_extract = pd.DataFrame({'g': np.array(record.annotations['abif_raw']['DATA9'])[index],
    'a': np.array(record.annotations['abif_raw']['DATA10'])[index],
    't': np.array(record.annotations['abif_raw']['DATA11'])[index],
    'c': np.array(record.annotations['abif_raw']['DATA12'])[index]})
    base_extract.to_csv(seq_file+'_extract.csv')
    handle.close()


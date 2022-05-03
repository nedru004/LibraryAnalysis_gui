#!/usr/bin/env python
import os
import numpy as np
import pandas as pd
import re
import itertools
import math
import multiprocessing as mp
import Bio
import os
import matplotlib.pyplot as plt
from plotnine import *
import seaborn as sns
import pymol

script_path = os.path.abspath(os.path.join(os.path.dirname(os.path.realpath(__file__)), '../bbmap/current/'))

codon = {'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L', 'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L', 'ATT': 'I',
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
    bbmerge_command += f' in1={in1} in2={in2} out1={out1} out2={out2} ecco mix'
    print(bbmerge_command)
    ret = os.system(bbmerge_command)
    assert ret == 0
    return out1, out2


def chunkit(align_files, file_size, chunk_size, other_args=[]):
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
            # Handle the case when a line is too long to fit the chunk size
            if chunk_start == chunk_end:
                chunk_end = get_next_line_position(chunk_end)
            # Save `process_chunk` arguments
            args = (align_files, chunk_start, chunk_end) + other_args
            chunk_args.append(args)
            # Move to the next chunk
            chunk_start = chunk_end
    return chunk_args


def david_call_parallel_variants(sam_file, wt, outfile, app, root):
    cpu_count = mp.cpu_count()
    file_size = os.path.getsize(sam_file)
    chunk_size = file_size // cpu_count
    chunk_args = chunkit(sam_file, file_size, chunk_size, other_args=[wt, outfile, app, root])
    # parallel
    with mp.Pool(cpu_count) as pool:
        answers = pool.starmap(david_call_variants, chunk_args)
    for answer in answers:
        correlations_wt = {k: correlations_wt.get(k, 0) + answer.get(k, 0) for k in set(correlations_wt) | set(answer)}
    wt_matrix = pd.DataFrame({'pair': list(correlations_wt.keys()), 'count': list(correlations_wt.values())})
    # # combine results
    # wt_count = None
    # for result in results:
    #     for xi, x in enumerate(result):
    #         wt_count[xi] += x
    #     totalcounts += result[0]
    # return totalcounts, wt_count


def david_call_variants(sam_file, wt, outfile, app, root):
    quality = int(app.quality.get())
    quality_nt = int(app.quality_nt.get())
    # Analyze Codons
    loop_read = open(sam_file, 'r')
    mutation_dict = {}
    variants = {'WT': 0}
    wt_count = [0]*int(len(wt)/3)
    file = open(outfile + '.csv', 'w')
    wt_file = open(outfile + '_wt.csv', 'w')
    variant_check = app.variant_check.get()
    correlation_check = app.correlations_check.get()
    read_count = 0
    percentage_reads = list(range(0, app.reads, int(app.reads/99)+(app.reads % 99 > 0)))
    while True:  # loop through each read from ngs
        try:
            tmp_r1 = next(loop_read).split('\t')
        except:
            break
        if '@' in tmp_r1[0]:
            continue
        # if read is paired to the next file
        if tmp_r1[6] == '=':
            tmp_r2 = next(loop_read).split('\t')
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
        tmpcounts = {tmp_r1[0]}  # initialize name
        tmp_wt_count = [0]*int(len(wt)/3+1)
        read_length = [100000, 0]
        for r3 in read_list:
            if (bin(int(r3[1]))[::-1][1] or bin(int(r3[1]))[::-1][3]) and int(r3[4]) > quality:  # if forward read is aligned and higher than quality threshold
                # set index
                ref = int(r3[3])-1  # subtract 1 for python numbering
                read_length[0] = min(read_length[0], ref)
                read = 0
                # cycle through each codon until end of read or end of reference sequence
                fragment_cigar = list(filter(None, re.split(r'(\d+)', r3[5])))
                for cigar in [fragment_cigar[ii:ii+2] for ii in range(0, len(fragment_cigar), 2)]:
                    if 'X' == cigar[1]:
                        # this is a substitution=
                        for i_sub in range(int(cigar[0])):  # not sure how this range works
                            codon_diff = (ref + i_sub) % 3  # find the correct frame
                            i_codon = r3[9][read + i_sub - codon_diff:read + i_sub - codon_diff + 3].upper()
                            # record codon if you can see the entire codon
                            if read + i_sub - codon_diff + 3 < len(r3[9]) and read + i_sub - codon_diff >= 0 and 'N' not in i_codon and '-' not in i_codon:
                                # only record if confident about entire codon
                                quality_codon = r3[10][read + i_sub - codon_diff:read + i_sub - codon_diff + 3]
                                if ord(quality_codon[0])-33 > quality_nt and ord(quality_codon[1])-33 > quality_nt and ord(quality_codon[2])-33>quality_nt:
                                    tmpcounts.add(str(int((ref+i_sub-codon_diff)/3+1))+'_' + i_codon)
                        ref += int(cigar[0])
                        read += int(cigar[0])
                    elif '=' == cigar[1]:
                        # this is a perfect match
                        codon_start = ref + ref % 3
                        codon_end = ref + int(cigar[0]) - (int(cigar[0]) % 3)
                        if codon_end - codon_start >= 3:  # if the perfect match spans an entire codon then count it
                            # make sure this doesnt exceed the length of the gene using min
                            wt_range = range(int(codon_start/3+1), min(int((codon_end-3)/3+1)+1, int(len(wt)/3)), 1)  # need to add 1 for range
                            for wt_i in wt_range:
                                tmp_wt_count[wt_i] = 1
                        ref += int(cigar[0])
                        read += int(cigar[0])
                    elif 'D' == cigar[1]:
                        # this is a deletion
                        tmpcounts.add(str(int(ref/3+1))+'_'+'del'+cigar[0])
                        # insert '-' for gap
                        r3[9] = r3[9][:read] + '-'*int(cigar[0]) + r3[9][read:]
                        r3[10] = r3[10][:read] + '-'*int(cigar[0]) + r3[9][read:]
                        read += int(cigar[0])
                        ref += int(cigar[0])
                    elif 'I' == cigar[1]:
                        # this is an insertion
                        insertion = r3[9][read:read+int(cigar[0])].upper()
                        tmpcounts.add(str(int(ref/3+1)) + '_' + 'ins' + insertion)
                        read += int(cigar[0])
                    elif 'S' == cigar[1]:
                        # soft clipping
                        read += int(cigar[0])
                    elif 'M' == cigar[1]:
                        # bad read (N) - ignore
                        ref += int(cigar[0])
                        read += int(cigar[0])
                    else:
                        raise TypeError('Cigar format not found')
                read_length[1] = max(read_length[1], ref)
        if len(tmpcounts) > 1:  # only record if mutations found
            # append data
            sorted_counts = sorted(tmpcounts, reverse=True)
            #totalcounts.append(sorted_counts)  # this is really slow. Instead write to file
            file.write(','.join(sorted_counts) + '\n')
            mutation = sorted_counts[1:]
            if app.count_indels.get() or (all(['ins' not in x for x in mutation]) and all(['del' not in x for x in mutation])):  # only recording if no indels were found
                for mut in sorted(mutation):
                    try:
                        mutation_dict[mut] += 1  # simple dictionary of mutations
                    except KeyError:
                        mutation_dict[mut] = 1
            if variant_check or correlation_check:  # and read_length[0] == 0 and read_length[1] >= len(wt):
                # convert to AA
                variant_mut = []
                for mut in mutation:
                    pos = int(mut.split('_')[0])-1
                    if mut.split('_')[1] in codon:
                        if codon[mut.split('_')[1]] != codon[wt.seq[pos*3:pos*3+3]]:
                            variant_mut.append(mut.split('_')[0] + '_' + codon[mut.split('_')[1]])
                    else:
                        variant_mut.append(mut)
                if variant_mut:
                    try:
                        variants['+'.join(variant_mut)] += 1
                    except KeyError:
                        variants['+'.join(variant_mut)] = 1
        else: #if len(tmpcounts) == 1 and variant_check and read_length[0] == 0 and read_length[1] >= len(wt):
            variants['WT'] += 1
        if sum(tmp_wt_count) > 1:  # and not indel:
            wt_positions = [xi for xi, x in enumerate(tmp_wt_count) if x]  # all positions wt codon was found
            for wt_idx in wt_positions:
                wt_count[wt_idx] += 1  # only add to total count after both reads have been processed
            wt_file.write(','.join(str(x) for x in wt_positions) + '\n')
        read_count += 1
        if read_count in percentage_reads:
            root.update_idletasks()
            app.progress['value'] += 1
            root.update()
    # Add wt to mutation_dict
    for i, wt_c in enumerate(wt_count[1:]):
        mutation_dict[str(i+1)+'_'+str(wt.seq[i*3:i*3+3])] = wt_c
    loop_read.close()
    file.close()
    wt_file.close()
    with open(outfile + '_mutlist.csv', 'w') as f:
        f.write('position,codon,AA,count\n')
        for mut in mutation_dict.keys():
            f.write(','.join(mut.split('_'))+',')
            try:
                f.write(codon[mut.split('_')[1]]+',')
            except KeyError:
                f.write(mut.split('_')[1]+',')
            f.write(str(mutation_dict[mut])+'\n')
    mutations = pd.DataFrame({'name': list(mutation_dict.keys()), 'count': list(mutation_dict.values())})
    mutations[['position', 'AA']] = mutations['name'].str.split('_', expand=True)
    mutations = mutations.replace({'AA': codon})
    mutations['position'] = pd.to_numeric(mutations['position'])
    matrix = mutations.pivot_table(index='AA', columns='position', values='count', aggfunc=sum)
    matrix.to_csv(outfile + '_matrix.csv')
    if variant_check or correlation_check:
        with open(outfile + '_variants.csv', 'w') as f:
            f.write('count,mutation\n')
            for mut in variants.keys():
                f.write(str(variants[mut]) + ',')
                f.write(mut + '\n')
    #os.remove(sam_file)
    #return mutation_table, correlation_matrix, correlation_matrix_wt


def david_paired_analysis(mut_files, wt_files, app, root):
    percentage_reads = list(range(0, app.reads, int(app.reads/99)))
    root.update_idletasks()
    app.progress['value'] = 0
    root.update()
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


def align_all_bbmap(sequencing_file, reference, sam_file, bbmap_script=f'java -cp {script_path} align2.BBMap',
                    max_gap=500, paired_sequencing_file=None):
    if paired_sequencing_file:
        command = f"{bbmap_script} ref={reference} in={sequencing_file} " \
                  f"in2={paired_sequencing_file} maxindel={max_gap} local outm={sam_file}"
    else:
        command = f"{bbmap_script} ref={reference} in={sequencing_file} " \
                  f"maxindel={max_gap} local outm={sam_file}"
    print(command)
    ret = os.system(command)
    assert ret == 0


def align_domain(app, root, sequencing_file, reference, domains_file, domain_out_file, bbmap_script=r'java -cp ../bbmap/current/ align2.BBMap', paired_sequencing_file=None):
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
    quality = int(app.quality.get())
    quality_nt = int(app.quality_nt.get())
    wt = Bio.SeqIO.read(reference, 'fasta').seq
    wt_count = [0] * int(len(wt) / 3 + 1)
    root.update_idletasks()
    app.progress['value'] = 0
    root.update()
    fraction = 100/len(domains)
    for domain in domains:
        #fixed_5p = domain.seq[:end5]
        #fixed_3p = domain.seq[end3:]
        domain_reference = os.path.join(os.path.dirname(domains_file), domain.id + '.fasta')
        with open(domain_reference, 'w') as outfile:
            outfile.write(f">{domain.id}\n{domain.seq}\n")
            outfile.close()
        # run alignment
        # if not os.path.exists(domain_out):
        #     if paired_sequencing_file:
        #         command = f"{bbmap_script} ref={domain_reference} in={sequencing_file} in2={paired_sequencing_file} outm={domain_out} local vslow minid=0 k=8 maxindel=10 startpad=8000 stoppad=8000 tipsearch=10 padding=10"
        #     else:
        #         command = f"{bbmap_script} ref={domain_reference} in={sequencing_file} outm={domain_out} local vslow minid=0 k=8 maxindel=10 startpad=8000 stoppad=8000 tipsearch=10 padding=10"
        #     print(command)
        #     ret = os.system(command)
        #     assert ret == 0
        # os.remove(domain_reference)
        print('##### Retaining reads that match insert #####')
        insert_filtered_fnames = os.path.join(os.path.dirname(domains_file), domain.id + '_insert_filtered.fastq')
        command = [f'java -cp {script_path} jgi.BBDuk',
                   'ref=%s' % domain_reference,
                   'in=%s' % sequencing_file,
                   'in2=%s' % paired_sequencing_file,
                   'outm=%s' % insert_filtered_fnames,
                   'mink=%d' % end5,
                   'hdist=0 mm=f']
        print('Running: \n\t%s' % ' '.join(command))
        ret = os.system(' '.join(command))
        assert ret == 0, 'bbduk.sh failed for %s!' % domain_reference

        print('##### Retaining reads also matching backbone #####')
        insert_bbone_filtered_fnames = os.path.join(os.path.dirname(domains_file), domain.id + '_insert_bbone_filtered.fastq')
        command = [f'java -cp {script_path} jgi.BBDuk',
                   'ref=%s' % reference,
                   'in=%s' % insert_filtered_fnames,
                   'outm=%s' % insert_bbone_filtered_fnames,
                   'k=%d' % 12,
                   'hdist=1 mm=f']
        print('Running: \n\t%s' % ' '.join(command))
        ret = os.system(' '.join(command))
        assert ret == 0, 'bbduk.sh failed for %s!' % reference

        print('##### Masking insert in reads #####')
        filtered_masked_fnames = os.path.join(os.path.dirname(domains_file), domain.id + '_masked.fastq')
        command = [f'java -cp {script_path} jgi.BBDuk',
                   'ref=%s' % domain_reference,
                   'in=%s' % insert_bbone_filtered_fnames,
                   'out=%s' % filtered_masked_fnames,
                   'mink=%d' % end5,
                   'kmask=N hdist=1 mm=f']
        ret = os.system(' '.join(command))
        assert ret == 0, 'bbduk.sh failed for %s!' % domain_reference

        # print('##### Trimming insert from reads #####')
        # trimmed_5p_fnames = os.path.join(os.path.dirname(domains_file), domain.id + '_trimmed.filtered')
        # GrepTrimMulti(pattern_5p, filtered_masked_fnames, trimmed_5p_fnames,
        #               trim_after=True, trim_match=True)

        print('##### Aligning to backbone #####')
        aligned_fnames = os.path.join(os.path.dirname(domains_file), domain.id + '_filtered_trimmed_5p_aligned.sam')
        command = [f'java -cp {script_path} align2.BBMap',
                   'ref=%s' % reference,
                   'in=%s' % filtered_masked_fnames,
                   'out=%s' % aligned_fnames,
                   'touppercase=t nodisk overwrite minratio=0.1']
        ret = os.system(' '.join(command))
        assert ret == 0, 'bbmap.sh failed for %s!' % reference

        # write to csv and create table
        read_length = [100000, 0]
        file_out = open(domain_reference.split('.fasta')[0] + '.csv', 'w')
        if os.path.exists(aligned_fnames):
            with open(aligned_fnames, 'r') as file:
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
                        if (bin(int(r3[1]))[::-1][1] or bin(int(r3[1]))[::-1][3]) and int(r3[4]) > quality:
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
                                    # make sure this doesnt exceed the length of the gene using min
                                    wt_range = range(int(codon_start / 3) + 1, min(int(codon_end/3) + 1, int(len(wt) / 3)), 1)  # need to add 1 for range
                                    for wt_i in wt_range:
                                        wt_count[wt_i] += 1
                    if len(tmpcounts) > 1:  # only record if mutations found
                        # append data
                        sorted_counts = sorted(tmpcounts, reverse=True)
                        # totalcounts.append(sorted_counts)  # this is really slow. Instead write to file
                        file_out.write(','.join(sorted_counts) + '\n')
            os.remove(aligned_fnames)
        else:
            print(aligned_fnames+' Failed Alignment!')
        os.remove(insert_filtered_fnames)
        os.remove(insert_bbone_filtered_fnames)
        os.remove(filtered_masked_fnames)
        os.remove(domain_reference)
        root.update_idletasks()
        app.progress['value'] += fraction
        root.update()
    for wt_pos, count in enumerate(wt_count):
        domain_counts.loc[wt_pos, 'wt'] = count
    # convert NaN to 0
    domain_counts = domain_counts.fillna(0)
    domain_counts.to_csv(domain_out_file.split('.sam')[0]+'_table.csv')
    file_out.close()
    return domain_counts


def map_structure(structure_file, matrix, min_q, max_q, name):
    # pymol.finish_launching()
    pymol.cmd.load(structure_file)
    structure_name = os.path.basename(structure_file).split('.')[0]
    pymol.stored.list = []
    pymol.cmd.iterate("(name ca)", "stored.list.append(resi)")
    matrix2 = matrix.iloc[:, 2:].reindex(columns=[int(x) for x in pymol.stored.list], fill_value=0).fillna(0)
    pymol.stored.all = [float(x) for x in list(matrix2.iloc[21, :])]  # average row
    pymol.cmd.alter(selection=structure_name, expression='b=0.0')
    pymol.cmd.alter(selection=structure_name + ' and n. CA', expression='b=stored.all.pop(0)')
    pymol.cmd.spectrum(expression='b', palette='red_white_green', selection=structure_name + ' and n. CA',
                       minimum=min_q, maximum=max_q)
    pymol.cmd.save(name + '_mapping.pse')

def calc_enrichment(base, selected, name, mincount, structure_file):
    base = base.groupby([x for x in base.columns if x not in ['count', 'codon']], as_index=False).agg('sum')
    selected = selected.groupby([x for x in selected.columns if x not in ['count', 'codon']], as_index=False).agg('sum')
    # combine into one dataframe
    df = pd.merge(base, selected, left_on=[x for x in base.columns if x != 'count'], right_on=[x for x in selected.columns if x != 'count'], suffixes=('_base', '_select'), how='outer')
    df = df.fillna(0)
    if 'position' in list(df.columns):
        # by AA position
        for i in df.position.unique():
            selection = df.loc[df.position == i, :]
            inp_count = sum(selection.loc[:, 'count_base']+1)
            sel_count = sum(selection.loc[:, 'count_select']+1)
            for select, row in selection.iterrows():
                inp_score = selection.loc[select, 'count_base'] + 1
                sel_score = selection.loc[select, 'count_select'] + 1
                score = math.log(sel_score / sel_count) - math.log(inp_score / inp_count)
                se = math.sqrt((1 / sel_score) + (1 / inp_score) + (1 / sel_count) + (1 / inp_count))
                df.loc[select, 'score'] = score
                df.loc[select, 'std_error'] = se
        df.to_csv(name, index=False)
        df['mutation'] = df['position'].astype(str) + '_' + df['AA']
        # create matrix
        df = df.loc[(df['count_base'] > mincount) | (df['count_select'] > mincount),]
        matrix = df.pivot_table(index='AA', columns='position', values='score', aggfunc=np.mean)
        matrix['Average'] = matrix.mean(axis=1)
        matrix = matrix[['Average'] + [c for c in matrix if c not in ['Average']]]
        matrix.loc['AveragePosition'] = matrix.mean(axis=0)
        matrix.to_csv(name + '_matrix.csv')
        # align to structure
        if structure_file:
            quartile = matrix.loc['AveragePosition'].quantile([0.10, 0.90])
            map_structure(structure_file, matrix, quartile[0.1], quartile[0.9], name)
    else:
        # by variant
        inp_count = sum(df.loc[:, 'count_base'] + 1)
        sel_count = sum(df.loc[:, 'count_select'] + 1)
        for select, row in df.iterrows():
            inp_score = df.loc[select, 'count_base'] + 1
            sel_score = df.loc[select, 'count_select'] + 1
            score = math.log(sel_score / sel_count) - math.log(inp_score / inp_count)
            se = math.sqrt((1 / sel_score) + (1 / inp_score) + (1 / sel_count) + (1 / inp_count))
            df.loc[select, 'score'] = score
            df.loc[select, 'std_error'] = se
        df.to_csv(name, index=False)
    # volcano plot
    p = ggplot(df) + aes(x='score', y='std_error') + geom_point() + theme_minimal() + geom_text(
        data=df.loc[df['score'] > 2]) + aes(x='score', y='std_error', label='mutation')
    print(p)
    print('Finished enrichment calculations')


def combine(dataset, combineBy, name, threshold, structure_file):
    df = dataset[0]
    if 'position' in list(df.columns):
        for data in dataset[1:]:
            df = pd.merge(df, data, on=['position', 'AA'])
        df['mutation'] = df['position'].astype(str) + '_' + df['AA']
    else:
        for data in dataset[1:]:
            df = pd.merge(df, data, on='mutation')
    # correlation plot first
    if combineBy:
        tmpdf = df.loc[(df[[x for x in df.columns if 'std_error' in x]] < threshold).apply(sum, axis=1) > len([x for x in df.columns if 'std_error' in x])/2]
    else:
        tmpdf = df.loc[(df[[x for x in df.columns if 'count' in x]] > threshold).apply(sum, axis=1) > len([x for x in df.columns if 'count' in x])/2]
    g = sns.PairGrid(tmpdf[[x for x in df.columns if 'score' in x]])
    g.map_diag(sns.histplot)
    g.map_offdiag(sns.regplot)
    plt.show()
    plt.pause(1)
    if not combineBy:
        # sum count_base
        df['count_combined_base'] = df[[x for x in df.columns if 'count_base' in x]].sum(axis=1)
        # sum count_select
        df['count_combined_select'] = df[[x for x in df.columns if 'count_select' in x]].sum(axis=1)
        df = df.fillna(0)
        if 'position' in list(df.columns):
            # by AA position
            for i in df.position.unique():
                selection = df.loc[df.position == i, :]
                inp_count = sum(selection.loc[:, 'count_combined_base'] + 1)
                sel_count = sum(selection.loc[:, 'count_combined_select'] + 1)
                for select, row in selection.iterrows():
                    inp_score = selection.loc[select, 'count_combined_base'] + 1
                    sel_score = selection.loc[select, 'count_combined_select'] + 1
                    score = math.log(sel_score / sel_count) - math.log(inp_score / inp_count)
                    se = math.sqrt((1 / sel_score) + (1 / inp_score) + (1 / sel_count) + (1 / inp_count))
                    df.loc[select, 'combined_score'] = score
                    df.loc[select, 'combined_std_error'] = se
        else:
            # by variant
            inp_count = sum(df.loc[:, 'count_combined_base'] + 1)
            sel_count = sum(df.loc[:, 'count_combined_select'] + 1)
            for select, row in df.iterrows():
                inp_score = df.loc[select, 'count_combined_base'] + 1
                sel_score = df.loc[select, 'count_combined_select'] + 1
                score = math.log(sel_score / sel_count) - math.log(inp_score / inp_count)
                se = math.sqrt((1 / sel_score) + (1 / inp_score) + (1 / sel_count) + (1 / inp_count))
                df.loc[select, 'combined_score'] = score
                df.loc[select, 'combined_std_error'] = se
        df1 = df.loc[(df['count_combined_base'] > threshold) | (df['count_combined_select'] > threshold), ]
        df['log_count'] = np.log10(df['count_combined_base'] + df['count_combined_select'])
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
                    df.loc[index, 'combined_score'] = sum(score * (var_betaML+[x**2 for x in se]) ** -1) / sum((var_betaML+[x**2 for x in se]) ** -1)
                    df.loc[index, 'combined_std_error'] = np.sqrt(1 / sum([x**-2 for x in se]))
                    df.loc[index, 'combined_eps'] = eps
                else:
                    df.loc[index, 'combined_score'] = np.nan
                    df.loc[index, 'combined_std_error'] = np.nan
                    df.loc[index, 'combined_eps'] = np.nan
        df1 = df.loc[df['combined_std_error'] < threshold, ]
        # volcano plot
        p = ggplot(df) + aes(x='combined_score', y='combined_std_error') + geom_point() + theme_minimal() + \
            geom_text(data=df.loc[df['combined_score'] > 2]) + aes(x='combined_score', y='combined_std_error',
                                                                   label='mutation') + scale_y_reverse()
        print(p)
    df.to_csv(name+'/combined_results.csv', index=False)
    # create matrix
    if 'position' in df1.columns:
        matrix = df1.pivot_table(index='AA', columns='position', values='combined_score', aggfunc=np.mean)
        missing_positions = list(set(range(1, max(df['position'])+1)).difference(matrix.columns))
        matrix[missing_positions] = np.nan
        matrix['Average'] = matrix.mean(axis=1, skipna=True)
        matrix = matrix[['Average'] + [c for c in matrix if c not in ['Average']]]
        matrix.loc['AveragePosition'] = matrix.mean(axis=0, skipna=True)
        matrix.to_csv(name+'/combined_matrix.csv')
        # align to structure
        if structure_file:
            quartile = matrix.loc['AveragePosition'].quantile([0.10, 0.90])
            map_structure(structure_file, matrix, quartile[0.1], quartile[0.9], name)
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
        out2 = paired_sequencing_file.split('.fastq.gz')[0] + '_corrected.fastq.gz'
        command = [f'java -cp {script_path} jgi.BBMerge',
                   'in=%s' % sequencing_file,
                   'in2=%s' % paired_sequencing_file,
                   'out1=%s' % out1,
                   'out2=%s' % out2,
                   "adapters=%s" % adapters]
        ret = os.system(' '.join(command))
        assert ret == 0, 'finding barcodes failed!'
    with open(adapters, 'r') as file:
        next(file)
        left_adapter = next(file).strip()
        next(file)
        right_adapter = next(file).strip()
    tmpbarcodes = os.path.join(os.path.dirname(sequencing_file), 'tmp_barcodes.fastq')
    command = [f'java -cp {script_path} jgi.BBDuk',
               'in=%s' % sequencing_file,
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
               'ktrim=r k=23 mink=11 hdist=1 tpe tbo']
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

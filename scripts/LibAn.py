# python
import AlignmentAnalyze
import time
import Bio.SeqIO
import Bio.SeqRecord
import argparse
import os


def _count_generator(reader):
    b = reader(1024 * 1024)
    while b:
        yield b
        b = reader(1024 * 1024)


def main(raw_args=None):
    parser = argparse.ArgumentParser(description='Analyze sequencing data for mutations')
    parser.add_argument('-wt', '--wtseq', help='FASTA file containing the wildtype sequence', required=True)
    parser.add_argument('-sam', '--sam', help='SAM file containing the aligned reads', required=False)
    parser.add_argument('-s', '--seq', help='FASTQ file containing the sequencing reads', default=None)
    parser.add_argument('-p', '--paired', help='FASTQ file containing the paired sequencing reads', default=False)
    parser.add_argument('-d', '--domains', help='FASTA file containing the domains of the wildtype sequence')
    parser.add_argument('-m', '--muts', help='File containing the mutations to be analyzed')
    parser.add_argument('-a', '--aamuts', help='File containing the amino acid mutations to be analyzed')
    parser.add_argument('-o', '--output', help='Output file directory and name')
    parser.add_argument('-pb', '--pacbio', help='Use long read sequencing alignment', action='store_true')
    parser.add_argument('-v', '--variant', help='Variant Analysis', action='store_true')
    parser.add_argument('-vfull', '--variantfull', help='Full length Variant Analysis', action='store_true')
    parser.add_argument('-c', '--correlation', help='Correlation Analysis', action='store_true')
    parser.add_argument('-i', '--indel', help='Analyze indel mutations', action='store_true')
    parser.add_argument('-minq', '--minq', help='Minimum quality of reads to analyze', default=15)
    parser.add_argument('-minb', '--minb', help='Minimum quality of bases to analyze', default=15)
    parser.add_argument('-par', '--parallel', help='Run analysis with number of cores. Default is 1', type=int, default=1)
    parser.add_argument('-f', '--force', help='Force overwrite of existing files', action='store_true')
    parser.add_argument('-n', '--nuc', help='Analyze nucleotide mutations', action='store_true')
    parser.add_argument('-aa', '--aa', help='Analyze amino acid mutations', action='store_true')
    parser.add_argument('-indel', '--max_indel', help='Maximum indel size to analyze', default=20)

    args = parser.parse_args(raw_args)

    # variables from arguments
    assert os.path.isfile(args.wtseq), f"given refrence/wildtype file name '{args.wtseq}' does not exist!"
    wt_seq = Bio.SeqIO.read(args.wtseq, "fasta").upper()
    # Sanity checks
    assert Bio.SeqIO.read(args.wtseq, "fasta").seq.translate(), f"given refrence/wildtype sequence file " \
                                                                f"'{args.wtseq}' is not a valid FASTA file " \
                                                                f"containing one unambiguous DNA sequence!"
    assert args.domains is None or os.path.isfile(args.domains), f"given domains file, '{args.domains}', does not exist"
    assert args.muts is None or os.path.isfile(args.muts), f"given domains file, '{args.muts}', does not exist"
    assert args.aamuts is None or os.path.isfile(args.aamuts), f"given domains file, '{args.aamuts}', does not exist"

    sequencing_file = args.seq
    paired_sequencing_file = args.paired
    # output files
    if args.seq:
        rootname = os.path.splitext(args.seq)[0].split('.fastq')[0]
    elif args.sam:
        rootname = os.path.splitext(args.sam)[0].split('.sam')[0]

    # create a log file containing inputs, number of reads, number of aligned reads, and number of aligned reads with at least one mutation
    log = open(rootname + '_log.txt', 'w')
    log.write(f"Sequencing File: {args.seq} \n")
    log.write(f"Paired Sequencing File: {args.paired} \n")

    # initialize some variables
    programstart = time.time()

    print('Reference sequence read\n')

    if args.muts:
        args.muts_list = []
        file = open(args.muts, 'r')
        while True:
            line = file.readline().strip()
            codon = file.readline().strip()
            if not codon: break
            args.muts_list.append(''.join([i for i in line if i.isdigit()]) + '_' + codon)
        file.close()
    if args.aamuts:
        args.aamuts_list = []
        file = open(args.aamuts, 'r')
        for line in file:
            args.aamuts_list.append(line.strip())
        file.close()



    ### Process Sequencing Records ###
    # merge paired reads
    if paired_sequencing_file and not os.path.exists(rootname+'_corrected.fastq.gz'):
        print(f'Merging paired reads. ({os.path.basename(rootname)})\n')
        sequencing_file = AlignmentAnalyze.correct_pairs(sequencing_file, paired_sequencing_file)
    else:
        print('Sequencing files already merged. Using existing corrected files')

    if args.domains:
        print(f'Finding Domains\n')
        AlignmentAnalyze.align_domain(args, sequencing_file, args.wtseq, args.domains, rootname)

    # align reads (using bbmap)
    if not os.path.exists(f'{rootname}.sam'):  # and args.aamuts:
        message = f'Aligning all sequences from {sequencing_file} to {wt_seq} using bbmap.'
        print(message)
        print('Aligning sequencing reads to reference.\n')
        if args.parallel == 1:
            par = None
        else:
            par = args.parallel
        if not args.pacbio:
            AlignmentAnalyze.align_all_bbmap(sequencing_file, args.wtseq, f'{rootname}.sam', par, max_gap=len(wt_seq))
        else:
            AlignmentAnalyze.align_long_read(sequencing_file, args.wtseq, f'{rootname}.sam', par)
    else:
        print('Sequencing files already aligned. Using existing sam file')

    # determine the number of reads
    with open(f'{rootname}.sam', 'rb') as fp:
        c_generator = _count_generator(fp.raw.read)
        # count each \n
        reads = sum(buffer.count(b'\n') for buffer in c_generator)
        print('Total lines:', reads + 1)
    if args.paired:
        reads = int((reads-3) / 2)
    print('Total number of reads to analyze: '+str(reads))
    args.reads = reads

    # add number of reads to log
    log.write(f"Total number of aligned reads: {str(reads)} \n")

    # call variants / find mutations
    print(f'Calling mutations/variants\n')

    if args.parallel > 1:
        # multi core
        AlignmentAnalyze.david_call_parallel_variants(f"{rootname}.sam", wt_seq, rootname, args)
    else:
        # single core
        output = AlignmentAnalyze.david_call_variants(f"{rootname}.sam", 0, os.path.getsize(f"{rootname}.sam"), wt_seq, rootname, args)
        AlignmentAnalyze.write_output(output[0], output[1], output[2], wt_seq, rootname)
    if args.correlation:
        print(f'Finding paired mutations\n')
        AlignmentAnalyze.david_paired_analysis(rootname + '_variants.csv', rootname + '_wt.csv', args, log)
    # os.system(f'java -cp ../bbmap/current/ var2.CallVariants2 in={rootname}.sam ref={args.wtseq} ploidy=1 out={rootname}.vcf 32bit')


    seq_analyze_time = time.time() - programstart
    time_per_seq = round(seq_analyze_time / reads * 1000, 3)
    print(f'Analyzed {reads} in {round(seq_analyze_time, 1)} seconds. {time_per_seq} ms per read.')
    args.reads = reads
    print('Finished Analysis!\n')
    log.close()

if __name__ == "__main__":
    main()


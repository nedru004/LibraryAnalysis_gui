# python
import AlignmentAnalyze
import time
import Bio
import Bio.Entrez
import Bio.Seq
import Bio.SeqIO
import Bio.SeqRecord
import tkinter as tk
from tkinter import filedialog
from tkinter.ttk import Progressbar
import pandas as pd
import os
import re
import sys

def run(sequencing_file=None, paired_sequencing_file=None):

    # variables from arguments
    assert os.path.isfile(app.wt_file), f"given refrence/wildtype file name '{app.wtseq}' does not exist!"
    wt_seq = Bio.SeqIO.read(app.wt_file, "fasta").upper()
    if sequencing_file==None or paired_sequencing_file==None:
        sequencing_file = app.seq_file
        paired_sequencing_file = app.paired_file
    # Sanity checks
    assert os.path.isfile(sequencing_file), f"given sequencing file, '{sequencing_file}', does not exist"
    old_stdout = sys.stdout
    log_file = open('test.fasta'.split('.fa')[0] + "message.log", "w")
    sys.stdout = log_file
    assert Bio.SeqIO.read(app.wt_file, "fasta").seq.translate(), f"given refrence/wildtype sequence file " \
                                                                f"'{app.wt_file}' is not a valid FASTA file " \
                                                                f"containing one unambiguous DNA sequence!"
    assert app.domains_file is None or os.path.isfile(app.domains_file), f"given domains file, '{app.domains_file}', does not exist"
    assert app.muts_file is None or os.path.isfile(app.muts_file), f"given domains file, '{app.muts_file}', does not exist"
    assert app.aamuts_file is None or os.path.isfile(app.aamuts_file), f"given domains file, '{app.aamuts_file}', does not exist"
    assert (app.muts_file or app.domains_file or app.aamuts_file or not app.deep_check.get()), \
        "Targeted analysis cannot be done if neither mutations file nor domains files are provided"
    assert app.control_correlations_file is None or os.path.isfile(app.control_correlations_file), \
        f"given refrence correlations file, '{app.control_correlations_file}', does not exist"

    # output files
    rootname = os.path.splitext(sequencing_file)[0].split('.fastq')[0]

    # initialize some variables
    programstart = time.time()

    message = 'Reference sequence read\n'
    app.output.insert('end', message)
    root.update()

    if app.muts_file:
        app.muts_list = []
        file = open(app.muts_file, 'r')
        for line in file:
            app.muts_list.append(line.strip())
        file.close()
    if app.aamuts_file:
        app.aamuts_list = []
        file = open(app.aamuts_file, 'r')
        for line in file:
            app.aamuts_list.append(line.strip())
        file.close()

    ### Process Sequencing Records ###
    # merge paired reads
    if paired_sequencing_file and not os.path.exists(rootname+'_corrected.fastq.gz'):
        app.output.insert('end', f'Merging paired reads. ({os.path.basename(rootname)})\n')
        root.update()
        sequencing_file = AlignmentAnalyze.correct_pairs(sequencing_file, paired_sequencing_file)
    else:
        print('Sequencing files already merged. Using existing corrected files')

    # align reads (using bbmap)
    if not os.path.exists(f'{rootname}.sam'):  # and app.aamuts_file:
        message = f'Aligning all sequences from {sequencing_file} to {wt_seq} using bbmap.'
        print(message)
        app.output.insert('end', 'Aligning sequencing reads to reference.\n')
        root.update()
        AlignmentAnalyze.align_all_bbmap(sequencing_file, app.wt_file, f'{rootname}.sam', max_gap=len(wt_seq))
    else:
        print('Sequencing files already aligned. Using existing sam file')

    # determine the number of reads
    lines = sum(1 for i in open(f'{rootname}.sam', 'rb'))
    reads = int((lines-3) / 2)
    print('Total number of reads to analyze: '+str(reads))
    app.reads = reads
    app.output.insert('end', f'Total number of reads to analyze: {str(reads)}\n')
    root.update()

    # call variants / find mutations
    app.output.insert('end', f'Calling mutations/variants\n')
    root.update()
    # David's work
    AlignmentAnalyze.david_call_variants(f"{rootname}.sam", wt_seq, rootname, app, root)
    if app.correlations_check.get():
        app.output.insert('end', f'Finding paired mutations\n')
        root.update()
        AlignmentAnalyze.david_paired_analysis(rootname + '_variants.csv', rootname + '_wt.csv', app, root)
    # os.system(f'java -cp ../bbmap/current/ var2.CallVariants2 in={rootname}.sam ref={app.wt_file} ploidy=1 out={rootname}.vcf 32bit')

    if app.domains_file:
        app.output.insert('end', f'Finding Domains\n')
        root.update()
        AlignmentAnalyze.align_domain(app, root, sequencing_file, app.wt_file, app.domains_file, rootname)

    seq_analyze_time = time.time() - programstart
    time_per_seq = round(seq_analyze_time / reads * 1000, 3)
    print(f'Analyzed {reads} in {round(seq_analyze_time, 1)} seconds. {time_per_seq} ms per read.')
    app.reads = reads
    app.output.insert('end', f'Analyzed {reads} in {round(seq_analyze_time, 1)} seconds. {time_per_seq} ms per read.\n')
    app.output.insert('end', 'Finished Analysis!\n')
    app.progress['value'] = 100
    root.update()
    sys.stdout = old_stdout
    log_file.close()

def enrichment():
    AlignmentAnalyze.calc_enrichment(app.base, app.select, app.selectcount.split('.fastq')[0] + '_results.csv', int(app.mincount.get()), app.pdb)


def run_batch():
    # find paired_files
    batch_folder = filedialog.askdirectory(title="Select folder containing fastq files")
    paired_files = {}
    for file in os.listdir(batch_folder):
        if '.fastq' in file and 'corrected' not in file:
            try:
                paired_files[re.sub(r'_R\d_', '', file)].append(os.path.join(batch_folder, file))
            except KeyError:
                paired_files[re.sub(r'_R\d_', '', file)] = [os.path.join(batch_folder, file)]
    for key, value in paired_files.items():
        assert len(value) < 3, f"Too many files linked to '{value[0]}'"
        run(sequencing_file=value[0], paired_sequencing_file=value[1])
    print('Finished Batch Analysis')


def combine():
    AlignmentAnalyze.combine(app.datasets, app.combineBy.get(), app.dirname, float(app.mincount.get()), app.pdb)


class Application(tk.Frame):
    def __init__(self, master=None):
        super().__init__(master)
        self.master = master
        self.grid(column=0)
        self.create_widgets()

    def create_widgets(self):
        self.winfo_toplevel().title("Library Alignment and Analysis")
        self.quiet_check = tk.IntVar()
        self.verbose_check = tk.IntVar()
        self.deep_check = tk.IntVar()
        self.correlations_check = tk.IntVar()
        self.parallel_check = tk.IntVar()
        self.bryan = tk.IntVar()
        self.domain_check = tk.IntVar()
        self.AA_check = tk.IntVar()
        self.variant_check = tk.IntVar()
        self.count_indels = tk.IntVar()

        tk.Label(self, text='NGS alignment', font=("Arial", 16, 'bold')).grid(column=0)
        self.seq = tk.Button(self, text='Input sequencing reads (fastq)', command=self.browse_seq)
        self.seq.grid(column=0)
        self.seq_file = None
        self.seq_label = tk.Label(self)
        self.seq_label.grid(column=0)

        self.paired = tk.Button(self, text='Input paired sequencing reads (fastq)', command=self.browse_paired)
        self.paired.grid(column=0)
        self.paired_file = None
        self.paired_label = tk.Label(self)
        self.paired_label.grid(column=0)

        self.wtseq = tk.Button(self, text='Input wildtype sequence (fasta)', command=self.browse_wt)
        self.wtseq.grid(column=0)
        self.wt_file = None

        tk.Label(self, text='Minimum Read Quality Score').grid(column=0)
        self.quality = tk.Entry(self, textvariable=tk.StringVar(self, '20'))
        self.quality.grid(column=0)

        tk.Label(self, text='Minimum Base Quality Score').grid(column=0)
        self.quality_nt = tk.Entry(self, textvariable=tk.StringVar(self, '25'))
        self.quality_nt.grid(column=0)

        self.variant = tk.Checkbutton(self, text='variant analysis', variable=self.variant_check)
        self.variant.select()
        self.variant.grid(column=0)

        self.correlation = tk.Checkbutton(self, text='Correlation analysis', variable=self.correlations_check)
        self.correlation.grid(column=0)

        self.indel = tk.Checkbutton(self, text='InDel analysis', variable=self.count_indels)
        self.indel.grid(column=0)

        tk.Label(self, text='Other').grid(column=0)
        self.mut = tk.Button(self, text='designed mutation NT list (csv)', command=self.browse_mut)
        self.mut.grid(column=0)
        self.muts_file = None

        self.mutaa = tk.Button(self, text='designed mutation AA list (csv)', command=self.browse_mutaa)
        self.mutaa.grid(column=0)
        self.aamuts_file = None

        self.control_correlations = tk.Button(self, text='Correlation file - baseline (csv)', command=self.browse_correlation)
        self.control_correlations.grid(column=0)
        self.control_correlations_file = None

        self.domains = tk.Button(self, text='designed insertion domains (fasta)', command=self.browse_domains)
        self.domains.grid(column=0)
        self.domains_file = None

        tk.Label(self, text='Run', font='bold').grid(column=0)
        self.run = tk.Button(self, text='Run analysis', command=run)
        self.run.grid(column=0)
        tk.Label(self, text='Run Batch', font='bold)').grid(column=0)
        self.run_batch = tk.Button(self, text='Run batch analysis', command=run_batch)
        self.run_batch.grid(column=0)

        tk.Label(self, text='Count Barcodes').grid(column=0)
        self.barcode = tk.Button(self, text='Select fastq for barcode count', command=self.browse_barcode)
        self.barcode.grid(column=0)

        self.output = tk.Text(self, height=10, width=60)
        self.output.grid(column=0)

        self.progress = Progressbar(self, orient='horizontal', length=500, mode='determinate')
        self.progress.grid(column=0)

        tk.Label(self, text='Calculate Enrichment (select counts)', font=("Arial", 16, 'bold')).grid(column=1, row=0)
        self.base_button = tk.Button(self, text='Input unselected counts', command=self.browse_base)
        self.base_button.grid(column=1, row=1)
        self.basecount = None
        self.w = tk.Label(self)
        self.w.grid(column=1, row=2)
        self.select_button = tk.Button(self, text='Input selected counts', command=self.browse_select)
        self.select_button.grid(column=1, row=3)
        self.selectcount = None
        self.w2 = tk.Label(self)
        self.w2.grid(column=1, row=4)
        tk.Label(self, text='Threshold (SE or count)').grid(column=1, row=5)
        self.mincount = tk.Entry(self, textvariable=tk.StringVar(self, '10'))
        self.mincount.grid(column=1, row=6)
        self.pdb_button = tk.Button(self, text='Select PDB file', command=self.browse_pdb)
        self.pdb_button.grid(column=1, row=7)
        self.pdb = ''
        tk.Label(self, text='Run').grid(column=1, row=8)
        self.analyze = tk.Button(self, text='Run analysis', command=enrichment)
        self.analyze.grid(column=1, row=9)

        tk.Label(self, text='Combine Datasets', font=("Arial", 16, 'bold')).grid(column=1, row=12)
        self.combineBy = tk.IntVar()
        self.r1 = tk.Radiobutton(self, text="Combine by Counts", variable=self.combineBy, value=0)
        self.r1.grid(column=1, row=13)
        self.r2 = tk.Radiobutton(self, text="Combine by Score", variable=self.combineBy, value=1)
        self.r2.grid(column=1, row=14)
        self.data = tk.Button(self, text='Input dataset', command=self.browse_data)
        self.data.grid(column=1, row=15)
        self.datasets = []
        self.combine = tk.Button(self, text='Combine Datasets', command=combine)
        self.combine.grid(column=1, row=16)

    def browse_seq(self):
        self.seq_file = filedialog.askopenfilename(title="Select a File")
        if self.seq_file != '':
            self.seq.config(bg='green', activebackground='green', relief=tk.SUNKEN)
            self.seq_label.config(text=os.path.basename(self.seq_file))
        else:
            self.seq_file = None

    def browse_paired(self):
        self.paired_file = filedialog.askopenfilename(title="Select a File")
        if self.paired_file != '':
            self.paired.config(bg='green', activebackground='green', relief=tk.SUNKEN)
            self.paired_label.config(text=os.path.basename(self.paired_file))
        else:
            self.paired_file = None

    def browse_wt(self):
        self.wt_file = filedialog.askopenfilename(title="Select a File")
        if self.wt_file != '':
            self.wtseq.config(bg='green', activebackground='green', relief=tk.SUNKEN)
        else:
            self.wt_file = None

    def browse_mut(self):
        self.muts_file = filedialog.askopenfilename(title="Select a File")
        if self.muts_file != '':
            self.mut.config(bg='green', activebackground='green', relief=tk.SUNKEN)
        else:
            self.muts_file = None

    def browse_mutaa(self):
        self.aamuts_file = filedialog.askopenfilename(title="Select a File")
        if self.aamuts_file != '':
            self.mutaa.config(bg='green', activebackground='green', relief=tk.SUNKEN)
        else:
            self.aamuts_file = None

    def browse_domains(self):
        self.domains_file = filedialog.askopenfilename(title="Select a File")
        if self.domains_file != '':
            self.domains.config(bg='green', activebackground='green', relief=tk.SUNKEN)
        else:
            self.domains_file = None

    def browse_correlation(self):
        self.control_correlations_file = filedialog.askopenfilename(title="Select a File")
        if self.control_correlations_file != '':
            self.control_correlations.config(bg='green', activebackground='green', relief=tk.SUNKEN)
        else:
            self.control_correlations_file = None

    def browse_base(self):
        self.basecount = filedialog.askopenfilename(title="Select a File")
        if self.basecount != '':
            self.base_button.config(bg='green', activebackground='green', relief=tk.SUNKEN)
            self.base = pd.read_csv(self.basecount)
            self.w.config(text=os.path.basename(self.basecount))
        else:
            self.basecount = None

    def browse_select(self):
        self.selectcount = filedialog.askopenfilename(title="Select a File")
        if self.selectcount != '':
            self.select_button.config(bg='green', activebackground='green', relief=tk.SUNKEN)
            self.select = pd.read_csv(self.selectcount)
            self.w2.config(text=os.path.basename(self.selectcount))
        else:
            self.selectcount = None

    def browse_pdb(self):
        self.pdb = filedialog.askopenfilename(title="Select PDB file")
        if self.pdb != '':
            self.pdb_button.config(bg='green', activebackground='green', relief=tk.SUNKEN)
        else:
            self.pdb = None

    def browse_data(self):
        tmpfile = filedialog.askopenfilenames(title="Select results files")
        if tmpfile != '':
            for file in tmpfile:
                self.datasets.append(pd.read_csv(file))
                app.output.insert('end', os.path.basename(file)+'\n')
                self.dirname = os.path.dirname(file)
            root.update()

    def browse_barcode(self):
        self.barcode_file = filedialog.askopenfilenames(title="Select fastq files (R1 & R2)")
        if len(self.barcode_file) > 1:
            sequencing_file = self.barcode_file[0]
            paired_sequencing_file = self.barcode_file[1]
            # app.output.insert('end', 'Merging paired reads.\n')
            # root.update()
            # sequencing_file, paired_sequencing_file = AlignmentAnalyze.correct_pairs(sequencing_file, paired_sequencing_file)
        else:
            sequencing_file = self.barcode_file[0]
            paired_sequencing_file = None
        adapters = filedialog.askopenfilename(title="Select adapter fasta file")
        AlignmentAnalyze.find_barcodes(sequencing_file, paired_sequencing_file, adapters)


if __name__ == "__main__":
    root = tk.Tk()
    root.geometry('900x850')
    app = Application(master=root)
    app.mainloop()

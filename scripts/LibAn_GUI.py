# python
import AlignmentAnalyze
import LibAn
import tkinter as tk
from tkinter import filedialog
from tkinter.ttk import Progressbar
import pandas as pd
import os
import re
import sys


def run_LibAn():
    arguments = ['-wt', app.wt_file, '-minq', int(app.quality.get()), '-minb', int(app.quality_nt.get()), '-par', app.parallel.get()]
    # run through the rest of the arguments
    if app.seq_file:
        arguments.extend(['-s', app.seq_file[0]])
    if app.sam_file:
        arguments.extend(['-sam', app.sam_file])
    if app.paired_file:
        arguments.extend(['-p', app.paired_file])
    if app.variant_check.get():
        arguments.append('-v')
    if app.variantfull_check.get():
        arguments.append('-vfull')
    if app.count_indels.get():
        arguments.append('-i')
    if app.pacbio.get():
        arguments.append('-pb')
    if app.correlations_check.get():
        arguments.append('-c')
    if app.muts_file:
        arguments.extend(['-m', app.muts_file])
    if app.aamuts_file:
        arguments.extend(['-a', app.aamuts_file])
    if app.domains_file:
        arguments.extend(['-d', app.domains_file])
    # convert list to string
    arguments = [str(item) for item in arguments]
    print('To run Analysis, run the following command:')
    print('python LibAn.py ' + ' '.join(arguments))
    LibAn.main(arguments)


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
        arguments = ['-s', value[0], '-p', value[1], '-wt', app.wt_file,
                     '-minq', int(app.quality.get()), '-minb', int(app.quality_nt.get()), '-par', app.parallel.get()]
        # run through the rest of the arguments
        if app.variant_check.get():
            arguments.append('-v')
        if app.variantfull_check.get():
            arguments.append('-vfull')
        if app.count_indels.get():
            arguments.append('-i')
        if app.pacbio.get():
            arguments.append('-pb')
        if app.correlations_check.get():
            arguments.append('-c')
        if app.muts_file:
            arguments.append(['-m', app.muts_file])
        if app.aamuts_file:
            arguments.append(['-a', app.aamuts_file])
        if app.domains_file:
            arguments.append(['-d', app.domains_file])
        # convert list to string
        arguments = [str(item) for item in arguments]
        print('To run Analysis, run the following command:')
        print('python LibAn.py ' + ' '.join(arguments))
        LibAn.main(arguments)
    print('Finished Batch Analysis')


def enrichment():
    AlignmentAnalyze.calc_enrichment(app, root, app.base, app.select, app.selectcount.split('.csv')[0] + '_results.csv', int(app.mincount.get()), app.pdb, app.muts_file)


def combine():
    AlignmentAnalyze.combine(app, root, app.datasets, app.combineBy.get(), app.dirname, float(app.mincount.get()), app.pdb)
    app.datasets = []
    # clear output window
    app.output.delete('1.0', tk.END)


def sanger_analysis():
    AlignmentAnalyze.analyze_sanger(filedialog.askopenfilename(title="Select ab1 file"))

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
        self.parallel = tk.IntVar()
        self.bryan = tk.IntVar()
        self.domain_check = tk.IntVar()
        self.AA_check = tk.IntVar()
        self.variant_check = tk.IntVar()
        self.variantfull_check = tk.IntVar()
        self.count_indels = tk.IntVar()
        self.pacbio = tk.IntVar()
        self.multiply = tk.IntVar()

        tk.Label(self, text='NGS alignment', font=("Arial", 16, 'bold')).grid(column=0)
        self.seq = tk.Button(self, text='Input sequencing file (or files) (fastq)', command=self.browse_seq)
        self.seq.grid(column=0)
        self.seq_file = None
        self.seq_label = tk.Label(self)
        self.seq_label.grid(column=0)

        self.paired = tk.Button(self, text='Input paired sequencing file (fastq)', command=self.browse_paired)
        self.paired.grid(column=0)
        self.paired_file = False
        self.paired_label = tk.Label(self)
        self.paired_label.grid(column=0)

        self.wtseq = tk.Button(self, text='Input wildtype sequence (fasta)', command=self.browse_wt)
        self.wtseq.grid(column=0)
        self.wt_file = None

        self.sam = tk.Button(self, text='Aligned File (SAM)', command=self.browse_sam)
        self.sam.grid(column=0)
        self.sam_file = None

        tk.Label(self, text='Minimum Read Quality Score').grid(column=0)
        self.quality = tk.Entry(self, textvariable=tk.StringVar(self, '15'))
        self.quality.grid(column=0)

        tk.Label(self, text='Minimum Base Quality Score').grid(column=0)
        self.quality_nt = tk.Entry(self, textvariable=tk.StringVar(self, '15'))
        self.quality_nt.grid(column=0)

        self.pacbio_check = tk.Checkbutton(self, text='PacBio analysis', variable=self.pacbio)
        self.pacbio_check.grid(column=0)

        self.variant = tk.Checkbutton(self, text='Variant analysis', variable=self.variant_check)
        self.variant.select()
        self.variant.grid(column=0)
        self.variantfull = tk.Checkbutton(self, text='     Full length variant analysis', variable=self.variantfull_check)
        self.variantfull.select()
        self.variantfull.grid(column=0)

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

        self.domains = tk.Button(self, text='designed insertion domains (fasta)', command=self.browse_domains)
        self.domains.grid(column=0)
        self.domains_file = None

        tk.Label(self, text='Run', font='bold').grid(column=0)
        tk.Label(self, text='Number of computing threads (parallel)').grid(column=0)
        self.parallel = tk.Entry(self, textvariable=tk.StringVar(self, '1'))
        self.parallel.grid(column=0)

        self.run = tk.Button(self, text='Run analysis', command=run_LibAn)
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
        self.time_button = tk.Button(self, text='Input time series', command=self.browse_time)
        self.time_button.grid(column=1, row=5)
        self.timecount = None
        self.output_time = tk.Text(self, height=5, width=40)
        self.output_time.grid(column=1, row=6)

        tk.Label(self, text='Threshold (SE or count)').grid(column=1, row=7)
        self.mincount = tk.Entry(self, textvariable=tk.StringVar(self, '10'))
        self.mincount.grid(column=1, row=8)
        self.pdb_button = tk.Button(self, text='Select PDB file', command=self.browse_pdb)
        self.pdb_button.grid(column=1, row=9)
        self.pdb = ''
        tk.Label(self, text='Run').grid(column=1, row=10)
        self.analyze = tk.Button(self, text='Run analysis', command=enrichment)
        self.analyze.grid(column=1, row=11)

        tk.Label(self, text='Combine Datasets', font=("Arial", 16, 'bold')).grid(column=1, row=14)
        self.combineBy = tk.IntVar()
        self.r1 = tk.Radiobutton(self, text="Combine by Counts", variable=self.combineBy, value=0)
        self.r1.grid(column=1, row=15)
        self.r2 = tk.Radiobutton(self, text="Combine by Score", variable=self.combineBy, value=1)
        self.r2.grid(column=1, row=16)
        self.multiply_check = tk.Checkbutton(self, text='Multiply Score by Count', variable=self.multiply)
        self.multiply_check.grid(column=1, row=17)
        self.multiply_check.select()
        self.data = tk.Button(self, text='Input dataset', command=self.browse_data)
        self.data.grid(column=1, row=18)
        self.datasets = []
        self.combine = tk.Button(self, text='Combine Datasets', command=combine)
        self.combine.grid(column=1, row=19)

        tk.Label(self, text='Other', font=("Arial", 12, 'bold')).grid(column=1, row=22)
        self.sanger_button = tk.Button(self, text="Sanger Analysis", command=sanger_analysis)
        self.sanger_button.grid(column=1, row=23)

    def browse_seq(self):
        self.seq_file = list(filedialog.askopenfilenames(title="Select a File"))
        if self.seq_file != '':
            self.seq.config(bg='green', activebackground='green', relief=tk.SUNKEN)
            self.seq_label.config(text=os.path.basename(self.seq_file[0]))
        else:
            self.seq_file = None

    def browse_paired(self):
        self.paired_file = filedialog.askopenfilename(title="Select a File")
        if self.paired_file != '':
            self.paired.config(bg='green', activebackground='green', relief=tk.SUNKEN)
            self.paired_label.config(text=os.path.basename(self.paired_file))
        else:
            self.paired_file = False

    def browse_wt(self):
        self.wt_file = filedialog.askopenfilename(title="Select a File")
        if self.wt_file != '':
            self.wtseq.config(bg='green', activebackground='green', relief=tk.SUNKEN)
        else:
            self.wt_file = None

    def browse_sam(self):
        self.sam_file = filedialog.askopenfilename(title="Select a File")
        if self.sam_file != '':
            self.sam.config(bg='green', activebackground='green', relief=tk.SUNKEN)
        else:
            self.sam_file = None

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

    def browse_time(self):
        # loop until cancel
        self.timecount = filedialog.askopenfilename(title="Select a count File with first timepoint")
        count = 0
        self.base = []
        self.select = None
        self.selectcount = self.timecount
        while self.timecount != '':
            self.time_button.config(bg='green', activebackground='green', relief=tk.SUNKEN)
            tmp_dataset = pd.read_csv(self.timecount)
            # loop through tmp_dataset columns but not position or AA
            if 'position' in list(tmp_dataset.columns):
                col = tmp_dataset.columns[3]
                newcol = str(count) + '_' + col
                tmp_dataset.rename(columns={col: newcol}, inplace=True)
            else:
                col = tmp_dataset.columns[2]
                newcol = str(count) + '_' + col
                tmp_dataset.rename(columns={col: newcol}, inplace=True)
            self.base.append(tmp_dataset)
            # write filename to output_time
            self.output_time.insert('end', str(count) + '  :  ' + os.path.basename(self.timecount)+'\n')
            self.timecount = filedialog.askopenfilename(title="Select a count File with next timepoint")
            count += 1

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
                tmp_dataset = pd.read_csv(file)
                datasetname = os.path.basename(file).split('_')[0]
                #loop through tmp_dataset columns but not position or AA
                if 'position' in list(tmp_dataset.columns):
                    for col in tmp_dataset.columns[2:]:
                        newcol = datasetname + '_' + col
                        tmp_dataset.rename(columns={col: newcol}, inplace=True)
                else:
                    for col in tmp_dataset.columns[1:]:
                        newcol = datasetname + '_' + col
                        tmp_dataset.rename(columns={col: newcol}, inplace=True)
                self.datasets.append(tmp_dataset)
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
            paired_sequencing_file = False
        adapters = filedialog.askopenfilename(title="Select adapter fasta file")
        AlignmentAnalyze.find_barcodes(sequencing_file, paired_sequencing_file, adapters)


if __name__ == "__main__":
    root = tk.Tk()
    root.geometry('900x850')
    app = Application(master=root)
    app.mainloop()

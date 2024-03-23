import os
import subprocess
from joblib import delayed, Parallel
from .validate import kwargs2list, combine_kwargs

def trim_files(input, output, paired_end: bool, custom_kwargs: dict) -> None:
    if paired_end:
        # Preparatory work.
        parent = os.path.abspath(os.path.join(str(output), os.pardir))
        os.makedirs(parent + "/paired", exist_ok = True)
        os.makedirs(parent + "/unpaired", exist_ok = True)
        argsList = []
        
        # Sort files into 1's and 2's.
        file1s = {}
        file2s = {}
        for file in os.listdir(input.fastq):
            if ".1.fastq" in file:
                substr_start = file.find(".1.fastq")
                file1s[file[0:(substr_start - 1)]] = file
            elif ".2.fastq" in file:
                substr_start = file.find(".2.fastq")
                file2s[file[0:(substr_start - 1)]] = file
            else:
                raise ValueError("File %s is not part of a pair." % file)
        
        # Now go through all the files we have.
        basenames = set(file1s.keys()) | set(file2s.keys())
        for basename in basenames:
            try:
                file1 = file1s[basename]
                file2 = file2s[basename]
            except ValueError:
                raise ValueError("There must be a file that doesn't have "
                                    "both paired-end reads (.1.fastq and .2.fastq). "
                                    "Check your files.")
            argsList.append(["trimmomatic", "PE"] + kwargs2list(custom_kwargs) +
                                [f"{input.fastq}/{file1}",
                                 f"{input.fastq}/{file2}",
                                 f"{parent}/paired/{basename}.paired.1.fastq.gz",
                                 f"{parent}/unpaired/{basename}.unpaired.1.fastq.gz",
                                 f"{parent}/paired/{basename}.paired.2.fastq.gz",
                                 f"{parent}/unpaired/{basename}.unpaired.2.fastq.gz",
                                 "SLIDINGWINDOW:51:20"])
        Parallel(n_jobs=6)(delayed(subprocess.run)(jobArgSet) for jobArgSet in argsList)
    else:
        os.makedirs(str(output))
        argsList = []
        for file in os.listdir(input.fastq):
            argsList.append(["trimmomatic", "SE"] + kwargs2list(custom_kwargs) +
                            [f"{input.fastq}/{file}",
                             f"{output}/{file}",
                             "SLIDINGWINDOW:51:20"])
        Parallel(n_jobs=6)(delayed(subprocess.run)(jobArgSet) for jobArgSet in argsList)

def hisat2_align(input, output, paired_end: bool, custom_kwargs) -> None:
    if paired_end:
        files_1 = ",".join([os.path.join(input.fastq, file) for file in os.listdir(input.fastq)
                           if ".1.fastq" in file])
        files_2 = ",".join([os.path.join(input.fastq, file) for file in os.listdir(input.fastq)
                           if ".2.fastq" in file])
        subprocess.run(["hisat2"] +
                        kwargs2list(combine_kwargs({"-x": f"{input.idx}/ref",
                                                    "-1": f"{files_1}",
                                                    "-2": f"{files_2}",
                                                    "-S": f"{output}"},
                                                    custom_kwargs)))
    else:
        files = ",".join([os.path.join(input.fastq, file) for file in os.listdir(input.fastq)
                          if "fastq" in file])
        subprocess.run(["hisat2"] +
                        kwargs2list(combine_kwargs({"-x": f"{input.idx}/ref",
                                                    "-U": f"{files}",
                                                    "-S": f"{output}"},
                                                    custom_kwargs)))
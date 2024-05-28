import os
import subprocess
from joblib import delayed, Parallel
from .validate import kwargs2list, combine_kwargs

# Input: has properties fastq
# Output: has properties results
def trim_files(input: str, output: str, paired_end: bool, trimmer: str, 
               custom_kwargs: dict, threads: int = 6) -> None:
    if paired_end:
        # Preparatory work.
        parent = os.path.abspath(os.path.join(output, os.pardir))
        os.makedirs(parent + "/paired", exist_ok = True)
        os.makedirs(parent + "/unpaired", exist_ok = True)
        argsList = []
        
        # Sort files into 1's and 2's.
        file1s = {}
        file2s = {}
        for file in os.listdir(input):
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
                                [f"{input}/{file1}",
                                 f"{input}/{file2}",
                                 f"{parent}/paired/{basename}.paired.1.fastq.gz",
                                 f"{parent}/unpaired/{basename}.unpaired.1.fastq.gz",
                                 f"{parent}/paired/{basename}.paired.2.fastq.gz",
                                 f"{parent}/unpaired/{basename}.unpaired.2.fastq.gz",
                                 trimmer])
            print(" ".join(argsList[-1]))
        Parallel(n_jobs=threads)(delayed(subprocess.run)(jobArgSet) for jobArgSet in argsList)
    else:
        os.makedirs(output)
        argsList = []
        for file in os.listdir(input):
            argsList.append(["trimmomatic", "SE"] + kwargs2list(custom_kwargs) +
                            [f"{input}/{file}",
                             f"{output}/{file}",
                             trimmer])
        Parallel(n_jobs=6)(delayed(subprocess.run)(jobArgSet) for jobArgSet in argsList)

# output: has properties results and summaryFile.
def hisat2_align(input: str, index: str, output: str, paired_end: bool, custom_kwargs) -> None:
    if paired_end:
        files_1 = ",".join([os.path.join(input, file) for file in os.listdir(input)
                           if ".1.fastq" in file])
        files_2 = ",".join([os.path.join(input, file) for file in os.listdir(input)
                           if ".2.fastq" in file])
        args = ["hisat2"] + kwargs2list(combine_kwargs(custom_kwargs,
                                                       {"-x": f"{index}/ref",
                                                        "-1": f"{files_1}",
                                                        "-2": f"{files_2}",
                                                        "-S": f"{output}"}))
    else:
        files = ",".join([os.path.join(input, file) for file in os.listdir(input)
                          if "fastq" in file])
        args = ["hisat2"] + kwargs2list(combine_kwargs(custom_kwargs,
                                                       {"-x": f"{index}/ref",
                                                        "-U": f"{files}",
                                                        "-S": f"{output}"}))
    print(" ".join(args))
    subprocess.run(args)
# VALIDATION FUNCTIONS.

import os
import validators
import warnings
import urllib
import requests
import re
from typing import Collection

# Find groups among your data.
def find_groups(dir: str) -> Collection[str]:
    if not os.path.exists(dir): raise FileNotFoundError("Supplied data directory "
        "'%s' not found. Please change your configuration file." % dir)
    return [name for name in os.listdir(dir) if os.path.isdir(
        os.path.join(dir, name))]
    
# Find the control group, which is defined as the filename starting with the
# phrase 'control' (case-insensitive).
# Note that controls don't have to be used, so only error if more than one.
def find_control_group(group_list: Collection[str]) -> str:
    control_group = [name for name in group_list if name.upper().startswith("CONTROL")]
    if len(control_group) == 0: return None
    elif len(control_group) == 1: return control_group[0]
    else: raise RuntimeError("There are multiple control groups. Please ensure "
                             "controls are manually specified or that only one "
                             "group starts with the phrase 'control'.")

# Validates that you do not overwrite builtin arguments and then combines them. Only
# works for keyword arguments.
def combine_kwargs(builtin_kwargs: dict, custom_kwargs: dict) -> dict:
    # Error checking: if dicts are None, then make some.
    if builtin_kwargs is None: builtin_kwargs = {}
    if custom_kwargs is None: custom_kwargs = {}

    # Overlap and return merged dicts.
    overlap = set(builtin_kwargs.keys()).intersection(custom_kwargs.keys())
    if len(overlap) > 0:
        raise ValueError("Warning: cannot override builtin keyword arguments: " + str(overlap))
    return builtin_kwargs | custom_kwargs

def kwargs2str(kwargs: dict) -> str:
    return " ".join([" ".join([str(k), str(v)]) for k, v in kwargs.items()])

def kwargs2list(kwargs: dict) -> list[str]:
    out = []
    [out.extend([str(k), str(v)]) for k, v in kwargs.items()]
    return out

# Cannot handle redirects.
def get_filename_from_url(url: str, is_ucsc = False):
    # Warnings if can't determine the filename.
    fn = None
    if not is_ucsc:
        response = requests.get(url)
        if response.status_code < 200 or response.status_code > 299:
            raise RuntimeError(f"URL {url} doesn't lead anywhere. Check your config file.")
        if "content-disposition" in response.headers:
            match = re.search("filename=([^;]+);?", response.headers["content-disposition"])
            if match is not None:
                fn = match.group(1)
    if fn is None or not fn:
        if not is_ucsc:
            warnings.warn(f"No valid CONTENT-DISPOSITION header for {url}. "
                        "Filename might be wrong/lack the proper extension.")
        fn = os.path.split(url.rstrip("/"))[-1]
    return fn

def get_filename(arg: str, is_ucsc = False):
    if not isinstance(arg, str):
        raise ValueError(f"Custom genome arguments must be strings. Supplied {arg}.")
    if os.path.exists(arg):
        return os.path.split(arg.rstrip("/"))[-1]
    elif validators.url(arg):
        # Element is URL. First try for content disposition then do last slash.
        return get_filename_from_url(arg, is_ucsc=is_ucsc)
    else:
        raise ValueError(f"Path {arg} doesn't exist on disk or the web. "
                         "Check your config file.")

def is_ucsc_genome_string(arg: str):
    return isinstance(arg, str)

def get_ucsc_genome_parts(arg: str):
    return arg.split(".")

# returns FASTA and GTF locations in order.
def process_genome_argument(arg):
    if arg is None or not arg:
        ValueError("Argument 'genome' must be defined in config. It cannot be blank.")

    is_ucsc = is_ucsc_genome_string(arg) # is the argument a ucsc name? (excludes explicit UCSC links)
    if is_ucsc: #ucsc genome. use where possible
        args = get_ucsc_genome_parts(arg)
        if len(args) == 1:
            args.append("ncbiRefSeq") # default naming
        genome, naming = args
        fasta = f"https://hgdownload.soe.ucsc.edu/goldenPath/{genome}/bigZips/{genome}.fa.gz"
        gtf = f"https://hgdownload.soe.ucsc.edu/goldenPath/{genome}/bigZips/genes/{genome}.{naming}.gtf.gz"
    elif isinstance(arg, dict):
        fasta, gtf = arg["fasta"], arg["gtf"] # These arguments are mandatory together
    else:
        raise RuntimeError(f"Invalid type {type(arg)} supplied to be processed.")
    
    locs = [fasta, gtf]
    names = [get_filename(x, is_ucsc=is_ucsc) for x in locs]
    return [{"fn": name,
             "fp": loc,
             "is_url": validators.url(loc),
             "rule_fn": re.search("(.+?)(\\.gz(ip)?)?$", name).group(1),
             "gzipped": bool(re.search("\\.gz(ip)?$", name))} for loc, name in zip(locs, names)]

def process_homer_genome_argument(homer_genome_arg, genome_fasta_info, genome_gtf_info):
    if is_ucsc_genome_string(homer_genome_arg):
        return get_ucsc_genome_parts(homer_genome_arg)[0], None, False
    # Downloading separate FASTA and GTFs for only HOMER is unsupported
    # elif isinstance(homer_genome_arg, dict):
    #     locs = [homer_genome_arg["fasta"], homer_genome_arg["gtf"]]
    #     names = [get_filename(x) for x in locs]
    #     return [None] + [{"fn": name,
    #                       "fp": loc,
    #                       "is_url": validators.url(loc),
    #                       "rule_fn": re.search("(.+?)(\\.gz(ip)?)?$", name).group(1),
    #                       "gzipped": bool(re.search("\\.gz(ip)?$", name))} for loc, name in zip(locs, names)]
    elif homer_genome_arg is not None:
        raise RuntimeError(f"Invalid type {type(homer_genome_arg)} supplied to be processed.")
    # None falls through to keep same genome and gtf fasta info
    return ("genome/" + genome_fasta_info["rule_fn"],
            "genome/annotations/" + genome_gtf_info["rule_fn"],
            True)

    
    
    
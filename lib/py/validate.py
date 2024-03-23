# VALIDATION FUNCTIONS.

import os
import warnings
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
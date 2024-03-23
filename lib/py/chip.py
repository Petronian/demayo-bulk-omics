import subprocess
from joblib import Parallel, delayed
from re import compile, search
from .validate import combine_kwargs, kwargs2list

def plot_heatmap(input, output, heatmapGroups, custom_kwargs):
    argsList = []
    for i, fn in enumerate(heatmapGroups):
        # "All" is first group in the list so this works by making i [1,...,N]
        sortArgs = {"-m": str(input),
                    "-o": str(output[i])}
        if fn != "All": sortArgs["--sortUsingSamples"] = str(i)
        argsList.append(["plotHeatmap"] + kwargs2list(combine_kwargs(sortArgs, custom_kwargs)))
    Parallel(n_jobs=6)(delayed(subprocess.run)(jobArgSet) for jobArgSet in argsList)

def find_intersections(input, output, custom_kwargs):
    argList = []
    for i, _ in enumerate(input):
        temp = [str(x) for x in input]
        group = temp.pop(i)
        kwargs = combine_kwargs({"-wo": "",
                                    "-filenames": "",
                                    "-a": group,
                                    "-b": ""},
                                custom_kwargs)
        del kwargs["-b"] # do this weird hack to ensure -b and args remains on the end
        argList.append((["bedtools", "intersect"] + kwargs2list(kwargs) + ["-b"] + temp, str(output[i])))
    def find_intersections(argTuple):
        # Write some comments at the top of the file, and then run the bedtools intersect function.
        args, output = argTuple
        file = open(output, "w")
        file.write("# Command: %s\n" % " ".join(args))
        file.write("# See https://bedtools.readthedocs.io/en/latest/content/tools/intersect.html for column information.\n")
        file.close()
        file = open(output, "a")
        subprocess.run(args, text=True, stdout=file)
        file.close()
    Parallel(n_jobs=6)(delayed(find_intersections)(argTuple) for argTuple in argList)

def find_all_intersections(input, output, custom_kwargs):
    temp = [str(x) for x in input]
    kwargs = combine_kwargs({"-a": temp[0],
                             "-b": ""},
                            custom_kwargs)
    del kwargs["-b"] # do this weird hack to ensure -b and args remains on the end
    argsList = ["bedtools", "intersect"] + kwargs2list(kwargs) + ["-b"] + temp[1:]
    file = open(str(output), "w")
    subprocess.run(argsList, text=True, stdout=file)
    file.close()

# assumes the only things that will match the regex are treatment and control.
def find_treatment_control_depths(file):
    with open(file) as fd:
        curr_line = fd.readline()
        info = {}
        info_pattern = compile("# tags after filtering in (\w+): (\d+)")
        while curr_line and not len(info) >= 2:
            curr_match = search(info_pattern)
            if curr_match:
                info[curr_match.group(1).lower()] = curr_match.group(2).lower()
        if len(info) < 2:
            raise RuntimeError("Unable to determine sequencing depths for treatment and control.")
    return (info["treatment"], info["control"])


def pairwise_differential_peakcall(groups_treatment_bdg, groups_control_bdg, info_files, groups, output_dir, custom_kwargs):
    group_bdg_info = list(zip(groups_treatment_bdg, groups_control_bdg))
    argsList = []
    for i, info in enumerate(group_bdg_info):
        curr_treat, curr_control = info
        for j in range(i + 1, len(group_bdg_info)):
            other_treat, other_control = group_bdg_info[j]
            depth_treat, depth_control = find_treatment_control_depths(info_files[i])[-1], find_treatment_control_depths(info_files[j])[-1]
            argsList.append(["macs3", "bdgdiff"], + kwargs2list(combine_kwargs({"--t1": str(curr_treat),
                                                                                "--t2": str(other_treat),
                                                                                "--c1": str(curr_control),
                                                                                "--c2": str(other_control),
                                                                                "--d1": str(depth_treat),
                                                                                "--d2": str(depth_control),
                                                                                "--outdir": str(output_dir),
                                                                                "--prefix": "-".join([str(groups[i]), str(groups[j])])},
                                                                                custom_kwargs)))
    Parallel(n_jobs=6)(delayed(subprocess.run)(jobArgSet) for jobArgSet in argsList)
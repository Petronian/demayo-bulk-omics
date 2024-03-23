import subprocess
from joblib import Parallel, delayed
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
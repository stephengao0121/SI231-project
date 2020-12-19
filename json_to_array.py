import json
import pandas as pd
import numpy as np


def j2a(mpath, fpath):
    features = pd.read_csv(fpath, sep='\t', header=None).values.tolist()
    markers_dict, cell_types, gene_types = {}, [], []
    f = open(mpath, 'r')
    data = json.load(f)
    for i in data["data"]:
        if i["alias"] not in cell_types:
            cell_types.append(i["alias"])
    for i in features:
        gene_types.append(i[1])
    # construct the dictionary of genes and expression.
    for i in data["data"]:
        if i["gene"].upper() in gene_types:
            if i["gene"].upper() not in markers_dict.keys():
                markers_dict[i["gene"].upper()] = [0 for k in range(len(cell_types))]
            markers_dict[i["gene"].upper()][cell_types.index(i["alias"])] = np.exp(i["avg_diff"])
    feature = []
    for i in features:
        feature.append(i[1])
    return markers_dict, cell_types, feature

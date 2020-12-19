import numpy as np
import scanpy
import scipy
import json_to_array as ja


def dim_reduction(a, marker, d):
    u, s, vh = np.linalg.svd(a)
    ud = u[:, 0: d]
    a_reduced = np.transpose(ud).dot(a)
    p_marker = np.transpose(ud).dot(marker)
    return a_reduced, p_marker


def main():
    marker_matrix, data_matrix = [], []
    # construct the marker matrix.
    marker_path = "4-markers\\hcl_top_markers_Adult-Brain.json"
    feature_path = "3-filteredFeatureMatrices\\V1_Human_Brain_Section_1_filtered_feature_bc_matrix\\features.tsv"
    marker_dict, cell_type, features = ja.j2a(marker_path, feature_path)
    for i in marker_dict.keys():
        if i in features:
            marker_matrix.append(marker_dict[i])
    marker_matrix = np.array(marker_matrix)

    # construct data matrix with needed rows.
    data = scanpy.read("3-filteredFeatureMatrices\\V1_Human_Brain_Section_1_filtered_feature_bc_matrix\\matrix.mtx")
    data = data.X.todense()
    for i in marker_dict.keys():
        data_matrix.append(data[features.index(i)].tolist()[0])
    data_matrix = np.array(data_matrix)

    # dim reduction and ls.
    ar, pm = dim_reduction(data_matrix, marker_matrix, 3)
    # f = np.linalg.lstsq(pm, ar, rcond=None)
    # f = np.linalg.pinv(pm).dot(ar)
    f = scipy.optimize.nnls(pm, ar)
    print(f)


if __name__ == '__main__':
    main()

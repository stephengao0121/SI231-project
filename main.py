import numpy as np
# import scanpy
from scipy.optimize import nnls
import json_to_array as ja
import argparse
import pandas as pd

def set_parser(parser):
    '''
    Set the parser

    Keyword Arguments:
        parser - argparse.ArugumentParser type, parser to be set
    '''
    parser.add_argument('-A', '--spatial_transcriptome_matrix', 
                        required=True, 
                        dest='Apath', 
                        help='The spatial transcriptome matrix, in csv format')

    parser.add_argument('-M', '--marker_matrix', 
                        required=True, 
                        dest='Mpath', 
                        help='The marker matrix, in csv format')

    parser.add_argument('-O', '--opath', 
                        required=True, 
                        dest='opath', 
                        help='output path, ending in .csv ')

    parser.add_argument('-d', '--number_of_dimensions', 
                        required=False, 
                        dest='d',
                        help='demensions to be retained in dimension reduction (30)', 
                        default=30)


def remove_specific_features(Amatrix, Bmatrix):
    '''
    Leaving only shared features in A and B matrix

    Keyword Arguments:
        Amatrix - matrix A, in csv format with rows representing features
        Bmatrix - matrix B, in csv format with rows representing features
    '''
    intersect_index = np.intersect1d(Amatrix.index, Bmatrix.index)
    Amatrix = Amatrix.loc[intersect_index]
    Bmatrix = Bmatrix.loc[intersect_index]

    return Amatrix, Bmatrix
    

def dim_reduction(a, marker, d):
    u, s, vh = np.linalg.svd(a)
    ud = u[:, 0: d]
    a_reduced = np.transpose(ud).dot(a)
    p_marker = np.transpose(ud).dot(marker)
    return a_reduced, p_marker


def main(args):
    # marker_matrix, data_matrix = [], []
    # construct the marker matrix.
    # marker_path = "4-markers\\hcl_top_markers_Adult-Brain.json"
    # feature_path = "3-filteredFeatureMatrices\\V1_Human_Brain_Section_1_filtered_feature_bc_matrix\\features.tsv"
    # marker_dict, cell_type, features = ja.j2a(marker_path, feature_path)
    # for i in marker_dict.keys():
    #     if i in features:
    #         marker_matrix.append(marker_dict[i])
    # marker_matrix = np.array(marker_matrix)

    # construct data matrix with needed rows.
    # data = scanpy.read("3-filteredFeatureMatrices\\V1_Human_Brain_Section_1_filtered_feature_bc_matrix\\matrix.mtx")
    # data = data.X.todense()
    # for i in marker_dict.keys():
    #     data_matrix.append(data[features.index(i)].tolist()[0])
    # data_matrix = np.array(data_matrix)

    Amatrix = pd.read_csv(args.Apath, index_col=0)
    Mmatrix = pd.read_csv(args.Mpath, index_col=0)

    Amatrix, Mmatrix = remove_specific_features(Amatrix, Mmatrix)

    # dim reduction and ls.
    ar, pm = dim_reduction(Amatrix.values, Mmatrix.values, args.d)
    Aprime = pd.DataFrame(data=ar, columns=Amatrix.columns)
    Mprime = pd.DataFrame(data=pm, columns=Mmatrix.columns)

    Fmatrix = pd.DataFrame(columns=Aprime.columns, index=Mprime.columns)
    # f = np.linalg.lstsq(pm, ar, rcond=None)
    # f = np.linalg.pinv(pm).dot(ar)
    for spot in Aprime.columns:
        Fmatrix.loc[:, spot], rnorm = nnls(Mprime.values, Aprime.loc[:,spot].values)

    Fmatrix = Fmatrix.agg(lambda c: c / np.sum(c))

    Fmatrix = Fmatrix.transpose()
    Fmatrix.to_csv(args.opath)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Fast and Robust Spatial Deconvolution')

    set_parser(parser)

    args = parser.parse_args()

    main(args)

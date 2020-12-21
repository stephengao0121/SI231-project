'''
Prepare the marker matrix from single cell data

author: panxq
email: panxq@shanghaitech.edu.cn
'''

import argparse
import numpy as np
import pandas as pd
import scanpy

def set_parser(parser):
    '''
    Set the parser

    Keyword Arguments:
        parser - argparse.ArugumentParser type, parser to be set
    '''
    parser.add_argument('-E', '--expression_matrix',
                        dest='expression_matrix', 
                        required=True, 
                        help='expression matrix to be extracted')

    parser.add_argument('-B', '--barcode', 
                        dest='barcodes', 
                        required=True, 
                        help='barcode for the expression matrix in tsv format')
    
    parser.add_argument('-F', '--feature', 
                        dest='features', 
                        required=True, 
                        help='Features for the expression matrix in tsv format')
    
    parser.add_argument('-C', '--clustering_result', 
                        dest='clusters', 
                        required=True, 
                        help='clustering result in csv format')
    
    parser.add_argument('-O', '-opath', 
                        dest='opath', 
                        required=True, 
                        help='output path')

def main(args):
    '''
    The main function

    Keyword Arguments:
        args - arguments to be used
    '''
    # read in data
    # barcodes and features read from tsv
    barcodes = pd.read_csv(args.barcodes, sep='\t', 
                           names=['barcode'])
    features = pd.read_csv(args.features, sep='\t', 
                           names=['ID', 'name', 'data_type'])
    clusters = pd.read_csv(args.clusters)
    cluster_names = clusters['Cluster'].unique()

    # expression matrix read from mtx files
    expression_matrix = scanpy.read(args.expression_matrix)
    expression_matrix = expression_matrix.X.todense()
    # expression matrix annotated by barcodes and features
    expression_matrix = pd.DataFrame(data=expression_matrix, 
                                     index=features['ID'].values, 
                                     columns=barcodes.values[:, 0])
    marker_matrix = pd.DataFrame(index=features['ID'].values)

    # computing mean over clusters to obtain markers
    for name in cluster_names:
        cluster_barcodes =  clusters.loc[clusters['Cluster'] == name, 'Barcode'].values
        marker_matrix.loc[:, name] = expression_matrix.loc[:, cluster_barcodes].agg(np.mean, axis=1).values

    marker_matrix.to_csv(args.opath)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Prepare the marker matrix from single cell data')

    set_parser(parser)

    args = parser.parse_args()

    main(args)

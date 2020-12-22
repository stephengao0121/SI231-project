import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import argparse
import matplotlib as mpl

def set_parser(parser):
    '''
    Set the parser

    Keyword Arguments:
        parser - argparse.ArugumentParser type, parser to be set
    '''
    parser.add_argument('-F', '--fraction_matrix', 
                        required=True, 
                        dest='Fpath', 
                        help='The fraction matrix path')

    parser.add_argument('-S', '--spatil_info', 
                        required=True, 
                        dest='SpatialInfo', 
                        help='csv file containing spatial information')

    parser.add_argument('-T', '--type', 
                        required=True, 
                        dest='type', 
                        help='cell type to be visualized')

    parser.add_argument('-O', '--opath', 
                        required=True, 
                        dest='opath', 
                        help='the output path for figure')

def main(args):
    '''
    main function
    '''
    spatial_info = pd.read_csv(args.SpatialInfo, 
                               names=['isUnder', 'x', 'y', 'figure_x', 'figure_y'], 
                               index_col=0)
    Fmatrix = pd.read_csv(args.Fpath, 
                          index_col=0)

    data = np.zeros([(spatial_info['x'].max()+1), (spatial_info['y'].max()+1)])

    for barcode, row in Fmatrix.iterrows():
        x = spatial_info.loc[barcode, 'x']
        y = spatial_info.loc[barcode, 'y']

        fraction = row[args.type]

        data[x,y] = fraction
    
    plt.style.use('dark_background')

    fig = plt.figure(figsize=(5,5))
    ax = fig.add_subplot(1,1,1)
    im = ax.imshow(data, cmap='gray', aspect=2)
    cbar = fig.colorbar(im)

    ax.set_xticks([])
    ax.set_yticks([])

    fig.savefig(args.opath)




if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='visualizer for F matrix')

    set_parser(parser)

    args = parser.parse_args()

    main(args)


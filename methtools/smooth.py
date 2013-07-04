#!/usr/bin/env python

import pandas as pd
import sys
import argparse

def smooth( options ):
    """
    
    """
    kwargs = {'center': options.center}
    if options.smooth_function == 'gaussian':
        kwargs.update({'std': options.gs})
    elif options.smooth_function == 'general_gaussian':
        kwargs = ({'power': options.ggp, 'width': options.ggw})
    elif options.smooth_function == 'kaiser':
        kwargs.update({'beta': options.kb})
    elif options.smooth_function == 'slepian':
        kwargs.update({'width': options.sw})

    print kwargs, options

    df = pd.read_csv( options.infile, sep='\t', index_col=None, names=['chr', 'start', 'end', 'name', 'value', 'strand'], dtype={'strand': object} )
    df['aggregate'] = df.groupby( ['chr'] ).apply(lambda x: pd.rolling_window(x['value'], options.window_length, options.smooth_function, **kwargs ) )

    df.to_csv(options.outfile, index=False, header=False, sep='\t', na_rep='0',
        cols=['chr', 'start', 'end', 'name', 'aggregate', 'strand'],
        float_format='%.2f'
        )


if __name__=='__main__':

    parser = argparse.ArgumentParser(
        description='Moving windows and smooting functions.')

    parser.add_argument("-i", "--infile", required=True, 
        type=argparse.FileType('r'),
        help="Path to the sample file.")

    parser.add_argument('-o', '--outfile', required=True, 
        type=argparse.FileType('w'),
        default=sys.stdout)

    parser.add_argument('-s', '--smooth-function', dest="smooth_function", default='boxcar', help='Smoothing function to use. [rolling mean]')
    parser.add_argument('-c', '--center', action='store_true', default=False, help='Center the smoothed values.')
    parser.add_argument('-w', '--window-len', type=int, dest="window_length", default=3, help='Window length if the smooth functions needs one. It is recommended to use a odd number.')

    parser.add_argument('--gs', type=float, help='Gaussian std.')
    parser.add_argument('--kb', type=float, help='Kaiser beta.')
    parser.add_argument('--ggp', type=int, help='General Gaussian power.')
    parser.add_argument('--ggw', type=int, help='General Gaussian width.')
    parser.add_argument('--sw', type=int, help='slepian width.')

    options = parser.parse_args()
    smooth( options )






#!/usr/bin/env python
# -*- coding: UTF-8 -*-

import argparse
import sys
from methtools.plot import main as plot_main
from methtools.destrand import main as destrand_main
from methtools.tiling import main as tiling_main
from methtools.filter import main as filter_main
from methtools.dmr import main as dmr_main
from methtools.calling import main as calling_main

def main():
    toolname = sys.argv[1].strip().lower()
    del sys.argv[1]

    if toolname == 'plot':
        plot_main()
    elif toolname == 'destrand':
        destrand_main()
    elif toolname == 'tiling':
        tiling_main()
    elif toolname == 'dmr':
        dmr_main()
    elif toolname == 'filter':
        filter_main()
    elif toolname == 'calling':
        calling_main()

if __name__ == '__main__':
    main()

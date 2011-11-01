#!/usr/bin/env python

import matplotlib
matplotlib.use('PDF')
from matplotlib import pyplot as plt


if __name__ == '__main__':
    import sys
    filenames = sys.argv[1:]
    from amsr_avhrr.validation import validate_all
    validate_all(filenames)
    #plt.show()

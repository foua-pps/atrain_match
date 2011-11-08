#!/usr/bin/env python

"""
Use this script to sum up results from a number of matched-values files.

"""
from validate_cph import validate_all


if __name__ == '__main__':
    from optparse import OptionParser
    parser = OptionParser(usage="Usage: %prog [options] FILE ...")
    parser.add_option('-v', '--verbose', action='store_true')
    parser.add_option('-l', '--land', action='store_true',
                      help="Process only land pixels")
    parser.add_option('-s', '--sea', action='store_true',
                      help="Process only sea pixels")
    opts, matched_values_files = parser.parse_args()
    
    if opts.land:
        if opts.sea:
            parser.error("land and sea options cannot be specified at the same "
                         "time")
        else:
            landsea = 'land'
    elif opts.sea:
        landsea = 'sea'
    else:
        landsea = None
    
    validate_all(matched_values_files, verbose=opts.verbose, landsea=landsea)

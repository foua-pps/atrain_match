'''
Created on Oct 27, 2010

@author: a001696
'''


def clean(crosses, file_finder=None):
    """Clean the SNO crosses in *crosses*, so that only those for which data is
    available are kept. If a file_finders object is provided, it will be used
    for searching for data files. Otherwise an HrptFileFinder will be used, with the
    time window in *crosses[0]* (or default HrptFileFinder time window if 
    crosses[0].time_window is None).
    
    Returns a dictionary with cross: data file pairs."""
    if file_finder is None:
        from file_finders import HrptFileFinder
        file_finder = HrptFileFinder()
        if crosses[0].time_window is not None:
            file_finder.set_time_window(crosses[0].time_window)
    
    cleaned = {}
    for cross in crosses:
        files_found = file_finder.find(cross.satellite1, cross.time1)
        if len(files_found) > 0:
            cleaned[cross] = files_found[0]
    
    return cleaned


def prepare_aapp(hrpt_filename, dir):
    """Create symbolic link to *hrpt_filename* in *aapp_working_dir/incoming*,
    in preparation for running AAPP to create Level 1B AVHRR files, suitable for
    PPS processing."""
    import os
    os.symlink(hrpt_filename, os.path.join(dir, os.path.basename(hrpt_filename)))


if __name__ == '__main__':
    from optparse import OptionParser
    
    # Set up and handle command line arguments
    parser = OptionParser()
    parser.set_usage("usage: %prog [options] FILE")
    parser.add_option('-l', '--link', type='string', metavar='DIR',
                      help="Add symbolic links to hrpt files in DIR, as a "
                      "preparation for running AAPP.")
    
    (options, args) = parser.parse_args()
    
    filename = args[0]
    
    from find_crosses import parse_crosses_file
    crosses = parse_crosses_file(filename)
    
    print("\n Satellites and Time Window:     %-12s%-12s%6.2f min\n" % \
          (crosses[0].satellite1, crosses[0].satellite2, crosses[0].time_window / 60.))
    for cross, hrpt_filename in sorted(clean(crosses).items()):
        print(cross.as_sno_line())
        try:
            prepare_aapp(hrpt_filename, options.link)
        except AttributeError:
            pass
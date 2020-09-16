#!/usr/bin/python
# -*- coding: utf-8 -*-
# Copyright (c) 2009-2019 atrain_match developers
#
# This file is part of atrain_match.
#
# atrain_match is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# atrain_match is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with atrain_match.  If not, see <http://www.gnu.org/licenses/>.
'''
Created on Oct 13, 2010

'''
from atrain_match.config import (RESOLUTION, _validation_results_dir, SURFACES)

import logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


def compile_stats(results_files, write=True, outfile_cfc="merged_sat_file_cfc", truth_sat='calipso'):
    """Run through all summary statistics."""

    # Always do CFC statistics, as it then that we read data
    # print("=========== Cloud fraction ============")
    from atrain_match.statistics import orrb_CFC_stat
    # read all results statistics only for cfc, resuse for cth, cty and cph
    cfc_stats = orrb_CFC_stat.CloudFractionStats(results_files=results_files, truth_sat=truth_sat)
    cfc_stats.write(outfile_cfc)

    if truth_sat not in ['amsr']:
        # note = "========== Cloud top height ==========="
        compiled_cth_file_name = outfile_cfc.replace('_cfc_', '_cth_')
        from atrain_match.statistics import orrb_CTH_stat
        cth_stats = orrb_CTH_stat.CloudTopStats(ac_data=cfc_stats.ac_data, truth_sat=truth_sat)
        cth_stats.write(compiled_cth_file_name)

    if truth_sat not in ['amsr']:
        # print("============= Cloud type ==============")
        compiled_cty_file_name = outfile_cfc.replace('_cfc_', '_cty_')
        from atrain_match.statistics import orrb_CTY_stat
        cty_stats = orrb_CTY_stat.CloudTypeStats(ac_data=cfc_stats.ac_data, truth_sat=truth_sat)
        cty_stats.write(compiled_cty_file_name)

    if truth_sat not in ['amsr']:
        # note = "========== Cloud phase ==========="
        compiled_cph_file_name = outfile_cfc.replace('_cfc_', '_cph_')
        from atrain_match.statistics import orrb_CPH_stat
        cth_stats = orrb_CPH_stat.CloudPhaseStats(ac_data=cfc_stats.ac_data, truth_sat=truth_sat)
        cth_stats.write(compiled_cph_file_name)

    # LWP only for AMSR-E currently
    if truth_sat in ['amsr']:
        # note = "========== Cloud lwp ==========="
        compiled_lwp_file_name = outfile_cfc.replace('_cfc_', '_lwp_')
        from atrain_match.statistics import orrb_LWP_stat
        cth_stats = orrb_LWP_stat.CloudLwpStats(ac_data=cfc_stats.ac_data, truth_sat=truth_sat)
        cth_stats.write(compiled_lwp_file_name)


if __name__ == '__main__':

    import os
    from glob import glob

    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--nodnt', '-d', const=True, nargs='?',
                        required=False,
                        help="Don't calculate statistics for all DNTs")
    parser.add_argument('--basic', '-b', const=True, nargs='?',
                        required=False,
                        help="Calculate the statistic for mode BASIC")
    parser.add_argument('--satz', '-z', const=True, nargs='?',
                        required=False,
                        help="Calculate the statistic for mode SATZ")
    parser.add_argument('--standard', '-t', const=True, nargs='?',
                        required=False,
                        help="Calculate the statistic for mode STANDARD")
    parser.add_argument('--odticfilter', '-o', const=True, nargs='?',
                        required=False,
                        help="Calculate the statistic for mode "
                        "OPTICAL_DEPTH_THIN_IS_CLEAR")
    parser.add_argument('--surface_new_way', '-n', const=True,
                        nargs='?', required=False,
                        help='calculate the statistic for the different '
                        'surfaces')
    parser.add_argument('--cotfilter', '-c', const=True,
                        nargs='?', required=False,
                        help='calculate the statistic for the optical '
                        'thickness filters')
    parser.add_argument('--satellites', '-sat', metavar='noaa18', type=str,
                        nargs='*', required=False,
                        help='List of satellite to combine, '
                        'overrides SATELLITES in atrain_match.cfg')
    parser.add_argument('--years', '-y', metavar='YYYY', type=str,
                        nargs='*', required=False,
                        help='List of year to combine, '
                        'overrides YEARS in atrain_match.cfg')
    parser.add_argument('--months', '-m', metavar='mm', type=str,
                        nargs='*', required=False,
                        help='List of year to combine, '
                        'overrides YEARS in atrain_match.cfg')      
    (options) = parser.parse_args()

    from atrain_match.utils.runutils import read_config_info
    AM_PATHS, SETTINGS = read_config_info()

    # Find all wanted modes (dnt)
    modes_list = []
    if options.nodnt:
        New_DNT_FLAG = ['']
    else:
        New_DNT_FLAG = ['', '_DAY', '_NIGHT', '_TWILIGHT']
    if options.basic:
        modes_list.append('BASIC')
    if options.standard:
        modes_list.append('STANDARD')
    if options.satz:
        for mode in ['SATZ_LOW', 'SATZ_70', 'SATZ_HIGH', 'SATZ_0_20', 'SATZ_20_40', 'SATZ_40_60', 'SATZ_60_80', 'SATZ_80_100',]:
            modes_list.append(mode)
    if options.odticfilter:  # I prefer this one! /KG
        print('Will calculate statistic for mode OPTICAL_DEPTH_THIN_IS_CLEAR')
        modes_list.append('OPTICAL_DEPTH_THIN_IS_CLEAR')
    if options.surface_new_way:
        print('Will calculate statistic for the different surfaces')
        for surface in SURFACES:
            modes_list.append(surface)
        # Add dnt flag to all modes so far
    modes_dnt_list = []
    for mode in modes_list:
        for dnt in New_DNT_FLAG:
            modes_dnt_list.append(mode + dnt)
    # Treat cotfilter separately as those directories have not dnt-flag at end
    if options.cotfilter:
        print('Will calculate statistic for mode COT-filter')
        for dnt in New_DNT_FLAG:
            for cot in SETTINGS["MIN_OPTICAL_DEPTH"]:
                # modes_list.append("OPTICAL_DEPTH-%0.2_%sf"(dnt, cot))# if like this
                modes_dnt_list.append("OPTICAL_DEPTH%s-%0.2f" % (dnt, cot))

    MY_YEARS = SETTINGS["YEARS"]
    MY_SATELLITES = SETTINGS["SATELLITES"]
    if options.satellites:
        MY_SATELLITES = options.satellites
    if options.years:
        MY_YEARS = options.years
    years_string = "_".join(MY_YEARS)    
    satellites_string = "_".join(MY_SATELLITES) 
    # get cases
    CASES = []
    month_list = ["*"]
    day_list = ["*"]
    if "MONTH" in SETTINGS.keys() and len(SETTINGS["MONTHS"]) > 0:
        month_list = ["{:02d}".format(int(ind)) for ind in SETTINGS["MONTHS"]]
    if options.months:
        month_list = ["{:02d}".format(int(ind)) for ind in options.months]
        
    if "DAY" in SETTINGS.keys() and len(SETTINGS["DAYS"]) > 0:
        day_list = ["{:02d}".format(int(ind)) for ind in SETTINGS["DAYS"]]
    for sat in MY_SATELLITES:
        for year in MY_YEARS:
            for month in month_list:
                for day in day_list:
                    CASES.append({'satname': sat,
                                  'month': month,
                                  'year': year,
                                  'day': day})

    if len(modes_dnt_list) == 0:
        logger.warning("No modes selected!")
        parser.print_help()
    # For each mode calcualte the statistics
    for process_mode_dnt in modes_dnt_list:
        print(RESOLUTION)
        # get result files for all cases
        for truth_sat in SETTINGS["COMPILE_STATISTICS_TRUTH"]:
            logger.info("PROCESS MODE %s, truth: %s", process_mode_dnt, truth_sat.upper())
            print("Gathering statistics from all validation results files in the "
                  "following directories:")
            results_files = []
            for case in CASES:
                indata_dir = AM_PATHS['result_dir'].format(
                    val_dir=_validation_results_dir,
                    satellite=case['satname'],
                    resolution=str(RESOLUTION),
                    month=case['month'],
                    year=case['year'],
                    day=case['day'],
                    mode=process_mode_dnt,
                    truth_sat=truth_sat,
                    min_opt_depth="")
                indata_dir = indata_dir.replace("_%H", "*")
                indata_file = AM_PATHS['result_file'].format(
                    resolution=str(RESOLUTION),
                    basename="*",
                    truth_sat=truth_sat)
                print("-> " + indata_dir)
                results_files.extend(glob("%s/*%skm*%s*.dat" % (indata_dir, RESOLUTION, truth_sat.lower())))
                results_files = list(set(results_files))
            if len(results_files) < 1:
                logger.info("PROCESS MODE %s have no results files", process_mode_dnt)
                continue
            # compile and write results
            compiled_dir = AM_PATHS['compiled_stats_dir'].format(
                val_dir=_validation_results_dir,
                satellite=satellites_string,
                resolution=str(RESOLUTION),
                month=case['month'],
                year=years_string,
                mode=process_mode_dnt,
                truth_sat=truth_sat)
            if not os.path.exists(compiled_dir):
                os.makedirs(compiled_dir)
            compiled_file_cfc = AM_PATHS['compiled_stats_filename'].format(
                satellite=satellites_string,
                resolution=str(RESOLUTION),
                month=case['month'],
                year=years_string,
                mode=process_mode_dnt,
                truth_sat=truth_sat,
                stat_type='cfc',
                min_opt_depth="")
            compiled_file_cfc = os.path.join(compiled_dir, compiled_file_cfc)
            compile_stats(results_files, outfile_cfc=compiled_file_cfc, truth_sat=truth_sat)

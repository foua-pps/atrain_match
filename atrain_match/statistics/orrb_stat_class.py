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

# Created on Oct 18, 2010

import numpy as np


class OrrbStats():
    """Abstract class for accumulating statistics from atrain_match."""

    def __init__(self, results_files=None, ac_data=None, truth_sat='calipso'):
        """Create an OrrbStats object with results from *results_files*."""
        self.results_files = results_files
        self.truth_sat = truth_sat
        self.ac_data = ac_data
        if self.results_files is not None:
            self.accumulate_data(results_files)
        if self.results_files is not None or ac_data is not None:
            self.do_stats()

    def read_one_file(self, datafile):
        # print(datafile)
        data_dict = {}
        current_datafile = open(datafile, "r")
        for line in current_datafile:
            if ":" not in line:
                continue
            line = line.replace('CALIOP', 'CALIPSO')
            # old before sept 2017 files have both cloudsat and calipso data in the same files
            # Note both CALIPSO and CALIOP where used
            if self.truth_sat.upper() not in line:
                continue
            what, data = line.rstrip().split(':')
            data = np.array([float(item) for item in data.lstrip().split(" ")])
            if what in data_dict.keys():
                print(what)
                raise KeyError("Key should not be already in list")
            data_dict[what] = data
        current_datafile.close()
        return data_dict

    def accumulate_data(self, results_files):
        print("reading data")
        acu = {}
        acu["scenes"] = len(results_files)
        # CFC DATA
        acu["n_clear_clear_cal"] = 0
        acu["n_clear_cloudy_cal"] = 0
        acu["n_cloudy_clear_cal"] = 0
        acu["n_cloudy_cloudy_cal"] = 0
        acu["n_clear_clear_cal_MODIS"] = 0
        acu["n_clear_cloudy_cal_MODIS"] = 0
        acu["n_cloudy_clear_cal_MODIS"] = 0
        acu["n_cloudy_cloudy_cal_MODIS"] = 0
        got_cloudsat_modis_flag = False
        # CTY

        acu["n_low_low"] = 0
        acu["n_low_medium"] = 0
        acu["n_low_high"] = 0
        acu["n_medium_low"] = 0
        acu["n_medium_medium"] = 0
        acu["n_medium_high"] = 0
        acu["n_high_low"] = 0
        acu["n_high_medium"] = 0
        acu["n_high_high"] = 0
        acu["n_cirrus_low"] = 0
        acu["n_cirrus_medium_tp"] = 0
        acu["n_cirrus_high_tp"] = 0
        acu["n_cirrus_medium_op"] = 0
        acu["n_cirrus_high_op"] = 0
        acu["n_clear_low"] = 0
        acu["n_clear_medium"] = 0
        acu["n_clear_high"] = 0
        acu["n_low_clear"] = 0
        acu["n_medium_clear"] = 0
        acu["n_high_clear"] = 0
        acu["n_frac_clear"] = 0
        acu["n_cirrus_clear"] = 0
        # CTH
        acu["cal_all_samples"] = {}
        acu["cal_low_samples"] = {}
        acu["cal_medium_samples"] = {}
        acu["cal_high_samples"] = {}
        acu["mean_error_cal_all_sum"] = {}
        acu["mean_error_cal_low_sum"] = {}
        acu["mean_error_cal_medium_sum"] = {}
        acu["mean_error_cal_high_sum"] = {}
        acu["rms_error_cal_all_sum"] = {}
        acu["rms_error_cal_low_sum"] = {}
        acu["rms_error_cal_medium_sum"] = {}
        acu["rms_error_cal_high_sum"] = {}
        acu["mae_error_cal_all_sum"] = {}
        acu["mae_error_cal_low_sum"] = {}
        acu["mae_error_cal_medium_sum"] = {}
        acu["mae_error_cal_high_sum"] = {}
        acu["n_over_250_cal_all"] = {}
        acu["n_over_250_cal_low"] = {}
        acu["n_over_250_cal_medium"] = {}
        acu["n_over_250_cal_high"] = {}
        acu["n_over_500_cal_all"] = {}
        acu["n_over_500_cal_low"] = {}
        acu["n_over_500_cal_medium"] = {}
        acu["n_over_500_cal_high"] = {}
        acu["n_over_1000_cal_all"] = {}
        acu["n_over_1000_cal_low"] = {}
        acu["n_over_1000_cal_medium"] = {}
        acu["n_over_1000_cal_high"] = {}
        acu["n_over_2500_cal_all"] = {}
        acu["n_over_2500_cal_low"] = {}
        acu["n_over_2500_cal_medium"] = {}
        acu["n_over_2500_cal_high"] = {}
        acu["n_missed_ctth_all"] = {}
        acu["n_missed_ctth_low"] = {}
        acu["n_missed_ctth_medium"] = {}
        acu["n_missed_ctth_high"] = {}
        acu["n_missed_cma_all"] = {}
        acu["n_missed_cma_low"] = {}
        acu["n_missed_cma_medium"] = {}
        acu["n_missed_cma_high"] = {}
        # CPH DATA
        acu["n_ice_ice_cal"] = 0
        acu["n_ice_water_cal"] = 0
        acu["n_water_ice_cal"] = 0
        acu["n_water_water_cal"] = 0
        # LWP
        acu["amsr_all_samples"] = 0
        acu["mean_error_amsr_all_sum"] = 0
        acu["rms_error_amsr_all_sum"] = 0

        cfc_stats_labels = ["CLOUD MASK %s-IMAGER TABLE" % (self.truth_sat.upper()),
                            "CLOUD MASK %s-PPS TABLE" % (self.truth_sat.upper())]
        cfcprob_stats_labels_clear = ["CLOUD MASK PROB %s-IMAGER TABLE CLEAR" % (self.truth_sat.upper())]
        cfcprob_stats_labels_cloudy = ["CLOUD MASK PROB %s-IMAGER TABLE CLOUDY" % (self.truth_sat.upper())]
        cfcprob_stats_labels_step = ["CLOUD MASK PROB %s-IMAGER TABLE STEP" % (self.truth_sat.upper())]

        cfc_stats_labels_modis = ["CLOUD MASK %s-MODIS TABLE" % (self.truth_sat.upper())]
        cty_stats_labels = ["CLOUD TYPE %s-IMAGER TABLE" % (self.truth_sat.upper()),
                            "CLOUD TYPE %s-PPS TABLE" % (self.truth_sat.upper())]
        cty_stats_labels_missed = [
            "CLOUD TYPE %s-IMAGER TABLE MISSED" % (self.truth_sat.upper()),
            "CLOUD TYPE %s-PPS TABLE MISSED" % (self.truth_sat.upper())]
        cph_stats_labels = ["CLOUD PHASE %s-IMAGER TABLE" % (self.truth_sat.upper())]
        lwp_stats_labels = ["CLOUD LWP %s-IMAGER TABLE" % (self.truth_sat.upper())]

        for datafile in self.results_files:
            data_dict = self.read_one_file(datafile)
            # Accumulate CALIOP/ISS/CLOUDSAT/AMSR-E statistics CFC
            for key in data_dict.keys():
                # If reprocessing old results files make sure to extract the right lines!
                # Ie do not use CLOUDSAT info when compiling stats for CALIPSO
                cal_data = data_dict[key]
                if key in cfc_stats_labels:
                    cal_data[cal_data < 0] = 0
                    acu["n_clear_clear_cal"] += cal_data[0]
                    acu["n_clear_cloudy_cal"] += cal_data[1]
                    acu["n_cloudy_clear_cal"] += cal_data[2]
                    acu["n_cloudy_cloudy_cal"] += cal_data[3]
                if key in cfcprob_stats_labels_step:

                    if "step_cmaprob" not in acu.keys():
                        acu["step_cmaprob"] = cal_data[0]
                    elif acu["step_cmaprob"] != cal_data[0]:
                        print("the same step in cma prob is not used for every file!")
                        raise ValueError
                if key in cfcprob_stats_labels_clear:
                    if "n_clear_cmaprob" not in acu.keys():
                        acu["n_clear_cmaprob"] = cal_data
                    else:
                        acu["n_clear_cmaprob"] = [acu["n_clear_cmaprob"][i] + cal_data[i] for i in range(len(cal_data))]
                if key in cfcprob_stats_labels_cloudy:
                    if "n_cloudy_cmaprob" not in acu.keys():
                        acu["n_cloudy_cmaprob"] = cal_data
                    else:
                        acu["n_cloudy_cmaprob"] = [acu["n_cloudy_cmaprob"][i] + cal_data[i]
                                                   for i in range(len(cal_data))]

                if key in cfc_stats_labels_modis:
                    got_cloudsat_modis_flag = True
                    modis_data = data_dict[key]
                    modis_data[modis_data < 0] = 0
                    acu["n_clear_clear_cal_MODIS"] += modis_data[0]
                    acu["n_clear_cloudy_cal_MODIS"] += modis_data[1]
                    acu["n_cloudy_clear_cal_MODIS"] += modis_data[2]
                    acu["n_cloudy_cloudy_cal_MODIS"] += modis_data[3]

            # Accumulate CALIOP statistics CPH
            for key in data_dict.keys():
                if key in cph_stats_labels:
                    cal_data = data_dict[key]
                    cal_data[cal_data < 0] = 0
                    acu["n_ice_ice_cal"] += cal_data[0]
                    acu["n_ice_water_cal"] += cal_data[1]
                    acu["n_water_ice_cal"] += cal_data[2]
                    acu["n_water_water_cal"] += cal_data[3]

            # Accumulate CALIOP statistics AMSRE
            for key in data_dict.keys():
                if key in lwp_stats_labels:
                    cal_data = data_dict[key]
                    acu["amsr_all_samples"] += cal_data[2]
                    acu["mean_error_amsr_all_sum"] += cal_data[2]*cal_data[0]
                    acu["rms_error_amsr_all_sum"] += cal_data[2]*cal_data[1]*cal_data[1]

            # Accumulate CALIOP/ISS/CLOUDSAT statistics CTY
            for key in data_dict.keys():
                if key in cty_stats_labels:
                    cal_data = data_dict[key]
                    cal_data[cal_data < 0] = 0
                    acu["n_low_low"] += cal_data[0]
                    acu["n_low_medium"] += cal_data[1]
                    acu["n_low_high"] += cal_data[2]
                    acu["n_medium_low"] += cal_data[3]
                    acu["n_medium_medium"] += cal_data[4]
                    acu["n_medium_high"] += cal_data[5]
                    acu["n_high_low"] += cal_data[6]
                    acu["n_high_medium"] += cal_data[7]
                    acu["n_high_high"] += cal_data[8]
                    acu["n_cirrus_low"] += cal_data[9]
                    acu["n_cirrus_medium_tp"] += cal_data[10]
                    acu["n_cirrus_high_tp"] += cal_data[11]
                    acu["n_cirrus_medium_op"] += cal_data[12]
                    acu["n_cirrus_high_op"] += cal_data[13]
                if key in cty_stats_labels_missed:
                    cal_data_missed = data_dict[key]
                    cal_data_missed[cal_data_missed < 0] = 0
                    acu["n_clear_low"] += cal_data_missed[0]
                    acu["n_clear_medium"] += cal_data_missed[1]
                    acu["n_clear_high"] += cal_data_missed[2]
                    acu["n_low_clear"] += cal_data_missed[3]
                    acu["n_medium_clear"] += cal_data_missed[4]
                    acu["n_high_clear"] += cal_data_missed[5]
                    acu["n_cirrus_clear"] += cal_data_missed[6]

            # Accumulate CALIOP/ISS/CLOUDSAT statistics CTH
            for key in data_dict.keys():
                if "HEIGHT" not in key:
                    continue
                type_of_clouds = key.split(" ")[-2]
                cloud_level = key.split(" ")[-1]
                data_one_cat = data_dict[key]
                if data_one_cat[3] < 0:
                    print("no pixels!")
                    continue
                if data_one_cat[3] <= 0:
                    data_one_cat[3] = 0

                if type_of_clouds not in acu["cal_{cloud_level}_samples".format(
                        cloud_level=cloud_level.lower())].keys():
                    acu["cal_{cloud_level}_samples".format(
                        cloud_level=cloud_level.lower())][type_of_clouds] = 0
                    acu["mean_error_cal_{cloud_level}_sum".format(
                        cloud_level=cloud_level.lower())][type_of_clouds] = 0
                    acu["rms_error_cal_{cloud_level}_sum".format(
                        cloud_level=cloud_level.lower())][type_of_clouds] = 0
                    acu["mae_error_cal_{cloud_level}_sum".format(
                        cloud_level=cloud_level.lower())][type_of_clouds] = 0
                    acu["n_missed_ctth_{cloud_level}".format(
                        cloud_level=cloud_level.lower())][type_of_clouds] = 0
                    acu["n_missed_cma_{cloud_level}".format(
                        cloud_level=cloud_level.lower())][type_of_clouds] = 0
                    acu["n_over_250_cal_{cloud_level}".format(
                        cloud_level=cloud_level.lower())][type_of_clouds] = 0
                    acu["n_over_500_cal_{cloud_level}".format(
                        cloud_level=cloud_level.lower())][type_of_clouds] = 0
                    acu["n_over_1000_cal_{cloud_level}".format(
                        cloud_level=cloud_level.lower())][type_of_clouds] = 0
                    acu["n_over_2500_cal_{cloud_level}".format(
                        cloud_level=cloud_level.lower())][type_of_clouds] = 0

                acu["n_missed_ctth_{cloud_level}".format(
                    cloud_level=cloud_level.lower())][type_of_clouds] += data_one_cat[5]
                acu["n_missed_cma_{cloud_level}".format(
                    cloud_level=cloud_level.lower())][type_of_clouds] += data_one_cat[4]
                acu["cal_{cloud_level}_samples".format(
                    cloud_level=cloud_level.lower())][type_of_clouds] += data_one_cat[3]
                acu["mean_error_cal_{cloud_level}_sum".format(
                    cloud_level=cloud_level.lower())][type_of_clouds] += data_one_cat[3]*data_one_cat[1]
                acu["rms_error_cal_{cloud_level}_sum".format(
                    cloud_level=cloud_level.lower())][type_of_clouds] += data_one_cat[3]*data_one_cat[2]*data_one_cat[2]
                acu["mae_error_cal_{cloud_level}_sum".format(
                    cloud_level=cloud_level.lower())][type_of_clouds] += data_one_cat[3]*data_one_cat[6]
                if len(data_one_cat) > 8:
                    acu["n_over_250_cal_{cloud_level}".format(
                        cloud_level=cloud_level.lower())][type_of_clouds] += data_one_cat[7]
                    acu["n_over_500_cal_{cloud_level}".format(
                        cloud_level=cloud_level.lower())][type_of_clouds] += data_one_cat[8]
                    acu["n_over_1000_cal_{cloud_level}".format(
                        cloud_level=cloud_level.lower())][type_of_clouds] += data_one_cat[9]
                    acu["n_over_2500_cal_{cloud_level}".format(
                        cloud_level=cloud_level.lower())][type_of_clouds] += data_one_cat[10]

        acu['got_cloudsat_modis_flag'] = got_cloudsat_modis_flag
        acu["Num"] = (acu["n_cloudy_cloudy_cal"] + acu["n_cloudy_clear_cal"] +
                      acu["n_clear_cloudy_cal"] + acu["n_clear_clear_cal"])
        self.ac_data = acu

    def do_stats(self):
        """
        Calculate all statistics and put results in instance attributes (?).

        """
        if not self.results_files and not self.ac_data:
            raise RuntimeError("No results files and no already loaded data")

    def printout(self):
        """Generate nice printout of the results."""
        raise NotImplementedError("The printout method should be implemented in"
                                  " a subclass of OrrbStat.")

    def write(self, filename, mode='w'):
        """Write printout to a file."""
        lines_to_write = self.printout()
        if len(lines_to_write) == 0:
            print("Not writing file %s" % (filename))
            print("No compiled results!")
        else:
            print("Writing file %s" % (filename))
            f = open(filename, mode)
            for l in lines_to_write:
                f.write(l + '\n')
            f.close()

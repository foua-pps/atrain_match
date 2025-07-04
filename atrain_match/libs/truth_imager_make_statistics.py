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
"""The main program for making statistics program."""


import atrain_match.config as config
# from atrain_match.config import INSTRUMENT
import os
import numpy as np
import logging

from atrain_match.utils.common import MatchupError

from atrain_match.libs.truth_imager_match import (get_matchups_from_data,
                                                  find_main_cloudproduct_file,
                                                  insert_info_in_filename_or_path,
                                                  add_additional_clousat_calipso_index_vars,
                                                  add_elevation_corrected_imager_ctth)

from atrain_match.truths.cloudsat import (add_validation_ctth_cloudsat,
                                          add_cloudsat_cloud_fraction)
from atrain_match.truths.calipso import (optical_depth_height_filtering,
                                         check_total_optical_depth_and_warn,
                                         add_validation_ctth_calipso,
                                         detection_height_filtering,
                                         set_thin_to_clear_filtering_1km)
from atrain_match.libs.truth_imager_statistics_lib import (calculate_statistics)
from atrain_match.plotting.trajectory_plotting import plot_satellite_trajectory
from atrain_match.plotting.along_track_plotting import (plot_cal_clsat_imager_time_diff,
                                                        plot_cal_clsat_geoprof_imager,
                                                        plot_cal_clsat_imager_satz,
                                                        plot_cal_clsat_cwc_imager)
from atrain_match.matchobject_io import (read_truth_imager_match_obj,
                                         CalipsoObject)

logger = logging.getLogger(__name__)


def add_validation_ctth(match_clsat, match_calipso):
    """ Add CPR (CloudSat) validation height if not in file already. """
    # Helper function to allow reprocessing of older matchup-files.
    if match_clsat is not None:
        if match_clsat.cloudsat.validation_height is None:
            match_clsat.cloudsat = add_validation_ctth_cloudsat(match_clsat.cloudsat)
        if match_clsat.cloudsat.cloud_fraction is None:
            match_clsat.cloudsat = add_cloudsat_cloud_fraction(match_clsat.cloudsat)
    if match_calipso is not None:
        if match_calipso.calipso.validation_height is None:
            match_calipso.calipso = add_validation_ctth_calipso(match_calipso.calipso)
    return match_clsat, match_calipso


def get_matchups(cross, AM_PATHS, SETTINGS, reprocess):
    """
    Retrieve matchups files. If *reprocess* is False, and if
    matchup files exist, get matchups directly from the processed files.
    Otherwise process matchup files first.
    """
    values = {}
    Obj_dict = {}
    out_dict = {}
    for truth in ['cloudsat', 'amsr', 'iss', 'synop', 'mora', 'calipso', 'dardar', 'earthcare']:
        Obj_dict[truth] = None
    try:
        values["satellite"] = cross.satellite1.lower()
    except AttributeError:
        raise ValueError('Need satellite1 and time (cross: %s)' % cross)

    if reprocess is False or SETTINGS['USE_EXISTING_RESHAPED_FILES']:
        for truth in ['cloudsat', 'amsr', 'iss', 'synop', 'mora', 'calipso', 'dardar', 'earthcare']:
            if not SETTINGS[truth.upper() + '_MATCHING']:
                logger.info(
                    "{truth} matching turned off {truth}_MATCHING]=False.".format(
                        truth=truth.upper()))
            else:
                values["atrain_sat"] = truth
                values["atrain_datatype"] = truth
                match_file, date_time = find_main_cloudproduct_file(
                    cross,
                    AM_PATHS['reshape_dir'],
                    AM_PATHS['reshape_file'],
                    values=values)
                if match_file is None:
                    logger.info(
                        "No processed {:s} match files found. ".format(truth) +
                        "Generating from source data if required.")
                    date_time = cross.time
                else:
                    Obj_dict[truth] = read_truth_imager_match_obj(match_file, truth=truth)
                    basename = '_'.join(os.path.basename(match_file).split('_')[1:5])

    if (all([obj_i is None for obj_i in Obj_dict.values()])):
        pass
    else:
        values['date_time'] = date_time
        values['year'] = date_time.year
        values['basename'] = basename
        values['month'] = "%02d" % (date_time.month)
        out_dict = {'basename': basename, 'values': values}

    # Redo matching if missing (we might want to remove this in future.)
    # And require that matchups are alredy done.
    # ================================================================
    if SETTINGS['USE_EXISTING_RESHAPED_FILES']:
        for truth in ['cloudsat', 'amsr', 'iss', 'synop',
                      'mora', 'calipso', 'dardar', 'earthcare']:
            if Obj_dict[truth] is None and SETTINGS[truth.upper()+'_REQUIRED']:
                raise MatchupError(
                    "Couldn't find calipso already processed matchup file, "
                    "USE_EXISTING_RESHAPED_FILES = True!")
    redo_matching = False
    for truth in ['cloudsat', 'amsr', 'iss', 'synop',
                  'mora', 'calipso', 'dardar', 'earthcare']:
        out_dict[truth] = Obj_dict[truth]
    if (all(obj_i is None for obj_i in Obj_dict.values())):
        out_dict = get_matchups_from_data(cross, AM_PATHS, SETTINGS)
    for truth in ['cloudsat', 'amsr', 'iss', 'synop',
                  'mora', 'calipso', 'dardar', 'earthcare']:
        if Obj_dict[truth] is None and SETTINGS[truth.upper()+'_REQUIRED']:
            redo_matching = True
    if redo_matching:
        out_dict = get_matchups_from_data(cross, AM_PATHS, SETTINGS)
    for truth in ['cloudsat', 'amsr', 'iss', 'synop',
                  'mora', 'calipso', 'dardar', 'earthcare']:
        if out_dict[truth] is None and SETTINGS[truth.upper()+'_REQUIRED']:
            raise MatchupError(
                "Couldn't find "
                "{truth} matchup and {truth}_REQUIRED is True!".format(
                    truth=truth))
    # ================================================================
    return out_dict


def plot_some_figures(match_clsat, match_calipso, values, basename, process_mode,
                      AM_PATHS, SETTINGS, am_obj=None, match_synop=None, mo_obj=None):
    """Plot the matchup track and trajectory etc. for one matchup file."""
    logger.info("Plotting")
    file_type = SETTINGS['PLOT_TYPES']

    plotpath = insert_info_in_filename_or_path(AM_PATHS['plot_dir'], values,
                                               datetime_obj=values['date_time'])
    if not os.path.exists(plotpath):
        os.makedirs(plotpath)

    # TRAJECTORY
    if match_calipso is not None and 1 == 2:
        imlon = match_calipso.imager.longitude.copy()
        imlat = match_calipso.imager.latitude.copy()
        trajectorypath = os.path.join(plotpath, "trajectory_plot")
        if not os.path.exists(trajectorypath):
            os.makedirs(trajectorypath)
        trajectoryname = os.path.join(trajectorypath,
                                      "%skm_%s_trajectory" % (int(config.RESOLUTION),
                                                              values['basename']))
        plot_satellite_trajectory(imlon,
                                  imlat,
                                  trajectoryname,
                                  config.AREA_CONFIG_FILE,
                                  file_type,
                                  **AM_PATHS)

    if (match_calipso is not None):
        # HEIGHT
        plot_cal_clsat_geoprof_imager(match_clsat,
                                      match_calipso,
                                      match_calipso.imager.imager_ctth_m_above_seasurface,
                                      plotpath,
                                      basename,
                                      process_mode,
                                      file_type,
                                      instrument=match_calipso.imager_instrument,
                                      MAXHEIGHT=SETTINGS["MAXHEIGHT"])
        # TIME DIFF SATZ
        plot_cal_clsat_imager_time_diff(match_clsat,
                                        match_calipso,
                                        plotpath, basename,
                                        config.RESOLUTION,
                                        instrument=match_calipso.imager_instrument)
        plot_cal_clsat_imager_satz(match_clsat,
                                   match_calipso,
                                   plotpath, basename,
                                   config.RESOLUTION, file_type,
                                   instrument=match_calipso.imager_instrument)

    if (match_clsat is not None and
            'rvod_liq_water_path' in match_clsat.cloudsat.all_arrays.keys()):

        elevation = np.where(np.less_equal(match_clsat.cloudsat.elevation, 0),
                             -9, match_clsat.cloudsat.elevation)
        data_ok = np.ones(match_clsat.cloudsat.elevation.shape, 'b')

        phase = 'LW'
        plot_cal_clsat_cwc_imager(match_clsat,
                                  elevation,
                                  data_ok,
                                  plotpath, basename,
                                  phase,
                                  instrument=match_clsat.imager_instrument)
        phase = 'IW'
        plot_cal_clsat_cwc_imager(match_clsat,
                                  elevation,
                                  data_ok,
                                  plotpath, basename, phase,
                                  instrument=match_clsat.imager_instrument)


def split_process_mode_and_dnt_part(process_mode_dnt):
    """Slit the process mode to two parts mode and day/night/twilight (DNT)."""
    mode_dnt = process_mode_dnt.split('_')
    if len(mode_dnt) == 1:
        process_mode = process_mode_dnt
        dnt_flag = None
    elif mode_dnt[-1] in ['DAY', 'NIGHT', 'TWILIGHT']:
        process_mode = '_'.join(mode_dnt[0:-1])
        dnt_flag = mode_dnt[-1]
    else:
        process_mode = process_mode_dnt
        dnt_flag = None
    return process_mode, dnt_flag


def process_one_mode(process_mode_dnt, match_calipso, match_clsat, iss_obj, am_obj, sy_obj, 
                     match_dardar, match_earthcare, min_optical_depth, values, AM_PATHS, SETTINGS, basename):
    """Make plots and statistics for one imager cloudproduct matchup (with one or several truths)."""

    # Get result filename
    process_mode, dnt_flag = split_process_mode_and_dnt_part(process_mode_dnt)
    min_depth_to_file_name = ""
    if process_mode == 'OPTICAL_DEPTH':
        min_depth_to_file_name = "-%.2f" % (min_optical_depth)
    values['mode'] = process_mode_dnt + min_depth_to_file_name
    result_path = insert_info_in_filename_or_path(AM_PATHS['result_dir'],
                                                  values,
                                                  datetime_obj=values['date_time'])
    if not os.path.exists(result_path):
        os.makedirs(result_path)
    result_file = AM_PATHS['result_file'].format(
        resolution=str(config.RESOLUTION),
        basename=values['basename'],
        truth_sat="xxx")
    statfilename = os.path.join(result_path, result_file)
    # Draw plots
    logger.debug("Plotting")
    if process_mode_dnt in SETTINGS['PLOT_MODES']:
        plot_some_figures(match_clsat, match_calipso, values, basename, process_mode,
                          AM_PATHS, SETTINGS, am_obj=am_obj)
    # Calculate Statistics
    logger.debug("Calculating statistics")
    calculate_statistics(process_mode, statfilename, match_calipso, match_clsat,
                         iss_obj, am_obj, sy_obj, match_earthcare, SETTINGS, dnt_flag)


def run(cross, run_modes, AM_PATHS, SETTINGS, reprocess=False):
    """The main work horse."""
    logger.info("Case: %s", str(cross))
    # sensor = INSTRUMENT.get(cross.satellite1.lower(), 'imager')

    if (not SETTINGS['USE_CMA_FOR_CFC_STATISTICS'] and
        not SETTINGS['USE_CT_FOR_CFC_STATISTICS'] and
            not SETTINGS['USE_CMAPROB_FOR_CFC_STATISTICS']):
        logger.error(
            "\n  --------------------------------- "
            "\n\tSet one of USE_*_FOR_CFC_STATISTICS=True in config.py!"
            "\n  --------------------------------- ")
        raise MatchupError("Configure problems, see messages above.")

    # Get the data that we need:
    matchup_results = get_matchups(cross, AM_PATHS, SETTINGS, reprocess)
    match_calipso = matchup_results['calipso']
    match_iss = matchup_results['iss']
    match_amsr = matchup_results['amsr']
    match_synop = matchup_results['synop']
    match_dardar = matchup_results['dardar']
    match_earthcare = matchup_results['earthcare']
    # mo_obj = matchup_results['mora']
    match_clsat = matchup_results['cloudsat']
    values = matchup_results['values']
    basename = matchup_results['basename']
    if match_calipso is not None and match_calipso.calipso.cloudsat_index is None:
        logger.info("Adding stuff missing in old reshaped files")
        match_clsat, match_calipso = add_additional_clousat_calipso_index_vars(match_clsat, match_calipso)
    logger.info("Adding validation height missing in old reshaped files")
    match_clsat, match_calipso = add_validation_ctth(match_clsat, match_calipso)
    # Calculate hight from sea surface
    match_clsat, match_calipso, match_iss, match_earthcare = add_elevation_corrected_imager_ctth(
        match_clsat, match_calipso, match_iss, match_earthcare, SETTINGS)
    calipso_original = CalipsoObject()
    # Save data orignal data that we might edit for some modes
    if match_calipso is not None:
        calipso_original.layer_top_altitude = match_calipso.calipso.layer_top_altitude.copy()
        calipso_original.layer_base_altitude = match_calipso.calipso.layer_base_altitude.copy()
        calipso_original.cloud_fraction = match_calipso.calipso.cloud_fraction.copy()
        calipso_original.feature_classification_flags = (
            match_calipso.calipso.feature_classification_flags.copy())
        calipso_original.validation_height = match_calipso.calipso.validation_height.copy()
        calipso_original.layer_top_pressure = match_calipso.calipso.layer_top_pressure.copy()
        calipso_original.layer_base_pressure = match_calipso.calipso.layer_base_pressure.copy()

    # For each mode, do the statistics:
    if (match_calipso is not None and
        SETTINGS['COMPILE_RESULTS_SEPARATELY_FOR_SINGLE_LAYERS_ETC'] and
        (SETTINGS['ALSO_USE_5KM_FILES'] or config.RESOLUTION == 5) and
            match_calipso.calipso.total_optical_depth_5km is None):
        logger.warning("\n\t Rematched_file is missing total_optical_depth_5km field"
                       "\n\t Consider reprocessing with: "
                       "\n\t COMPILE_RESULTS_SEPARATELY_FOR_SINGLE_LAYERS_ETC=True"
                       "\n\t ALSO_USE_5KM_FILES=True or RESOLUTION==5")

    for process_mode_dnt in run_modes:
        logger.info("Process mode: %s", process_mode_dnt)
        optical_depths = [None]         # Update this if you always want to do filtering!/Nina
        if process_mode_dnt in ["OPTICAL_DEPTH", "OPTICAL_DEPTH_DAY",
                                "OPTICAL_DEPTH_NIGHT", "OPTICAL_DEPTH_TWILIGHT"]:
            optical_depths = SETTINGS['MIN_OPTICAL_DEPTH']

        # split process_mode_dnt into two parts. One with process_mode and one dnt_flag
        process_mode, dnt_flag = split_process_mode_and_dnt_part(process_mode_dnt)
        for min_optical_depth in optical_depths:
            # For some modes these are updated, so reset calipso data to original
            if match_calipso is not None:
                # -------------------------------------------------------------
                match_calipso.calipso.layer_top_altitude = calipso_original.layer_top_altitude.copy()
                match_calipso.calipso.layer_base_altitude = calipso_original.layer_base_altitude.copy()
                match_calipso.calipso.cloud_fraction = calipso_original.cloud_fraction.copy()
                match_calipso.calipso.feature_classification_flags = calipso_original.feature_classification_flags.copy()
                match_calipso.calipso.validation_height = calipso_original.validation_height.copy()
                match_calipso.calipso.layer_top_pressure = calipso_original.layer_top_pressure.copy()
                match_calipso.calipso.layer_base_pressure = calipso_original.layer_base_pressure.copy()
                # -------------------------------------------------------------
            # If mode = OPTICAL_DEPTH -> Change cloud -top and -base profile
            if match_calipso is not None and process_mode == 'OPTICAL_DEPTH':
                use_old_method = SETTINGS['KG_OLD_METHOD_CLOUD_CENTER_AS_HEIGHT']
                retv = optical_depth_height_filtering(
                    match_calipso.calipso,
                    min_optical_depth,
                    use_old_method=use_old_method,
                    limit_ctop=SETTINGS['OPTICAL_LIMIT_CLOUD_TOP'])
                match_calipso.calipso.layer_top_altitude = retv[0]
                match_calipso.calipso.layer_base_altitude = retv[1]
                match_calipso.calipso.cloud_fraction = retv[2]
                match_calipso.calipso.feature_classification_flags = retv[3]
                match_calipso.calipso.validation_height = retv[4]
                match_calipso.calipso.layer_top_pressure = retv[5]
                match_calipso.calipso.layer_base_pressure = retv[6]
            if match_calipso is not None:
                check_total_optical_depth_and_warn(match_calipso)
                if 'STANDARD' in process_mode:
                    match_calipso.calipso.validation_height = detection_height_filtering(match_calipso)
            if match_calipso is not None and process_mode == 'OPTICAL_DEPTH_THIN_IS_CLEAR':
                logger.info("Setting thin clouds to clear, "
                            "using 5km data in mode OPTICAL_DEPTH_THIN_IS_CLEAR")
                retv = set_thin_to_clear_filtering_1km(match_calipso, SETTINGS)
                match_calipso.calipso.cloud_fraction = retv[0]
                match_calipso.calipso.validation_height = retv[1]
            # Time to process results files for one mode:
            process_one_mode(process_mode_dnt,
                             match_calipso, match_clsat, match_iss, match_amsr, match_synop, match_dardar,
                             match_earthcare,
                             min_optical_depth, values,
                             AM_PATHS, SETTINGS, basename)
    # We are done, free some memory:
    match_calipso = None
    match_clsat = None
    match_iss = None
    match_amsr = None
    match_dardar = None

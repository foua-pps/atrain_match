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
"""Read all matched data and make some plotting
"""
import os
import re
from glob import glob
import numpy as np

from utils.get_flag_info import get_calipso_clouds_of_type_i_feature_classification_flags_one_layer
import matplotlib.pyplot as plt
import matplotlib
from utils.get_flag_info import (get_calipso_low_medium_high_classification,
                           get_inversion_info_pps2014,
                           get_calipso_clouds_of_type_i,
                           get_calipso_high_clouds,
                           get_calipso_medium_clouds,
                           get_calipso_low_clouds)

from scipy.stats import kurtosis, skewtest, skew, mode, kurtosis

from utils.stat_util import (my_hist,
                       my_iqr,
                       my_rms,
                       my_mae,
                       half_sample_mode,
                       half_sample_mode,
                       my_pe250m,
                       my_pe500m,
                       my_pe1000m,
                       my_pe2000m,
                       my_pe2500m,
                       my_pe5000m)

matplotlib.rcParams.update(matplotlib.rcParamsDefault)
matplotlib.rcParams.update({'font.size': 20})
# plt.rc('text', usetex=True)
# plt.rc('font', family='serif')
# matplotlib.use('ps')
from matplotlib import rc
from my_dir import ADIR
# rc('text', usetex=True)
# rc('text.latex', preamble='\usepackage{color}')
out_filename = ADIR + "/Documents/A_PPS_v2017/Validation_2018/results_ctth_polarn.txt"
out_file_h = open(out_filename, 'a')


def my_make_plot3(y, x, x2, mhl, use):
    fig = plt.figure(figsize=(15, 11))
    ax = fig.add_subplot(111)
    use_k = np.logical_and(mhl["all_clouds_tp_thin_1layers"], use)
    abias1 = np.abs(y[use_k] - x[use_k])
    sort_ind = np.argsort(x[use_k])
    plt.plot(x[use_k][sort_ind], 'g.')
    plt.plot(y[use_k][sort_ind], 'b.', alpha=0.2)
    # plt.show()


def my_make_plot2(y, x, x2, mhl, use):
    fig = plt.figure(figsize=(15, 11))
    ax = fig.add_subplot(321)
    use_k = np.logical_and(use, x2 > 0)
    # use_k = np.logical_and(use_k, x-x2>3000)

    print(min(y[use_k]), len(use[use_k]))
    abias1 = np.abs(y[use_k] - x[use_k])
    abias2 = np.abs(y[use_k] - x2[use_k])
    dist = 0.5*np.abs(x2[use_k] - x[use_k])
    dist25 = 0.25*np.abs(x2[use_k] - x[use_k])
    closer_to_top = np.logical_and(abias1 <= abias2, np.logical_and(y[use_k] <= x[use_k], y[use_k] >= x2[use_k]))
    closer_to_2 = np.logical_and(abias1 > abias2, np.logical_and(y[use_k] <= x[use_k], y[use_k] >= x2[use_k]))
    in_between = np.logical_and(y[use_k] < (x[use_k]), y[use_k] > (x2[use_k]))
    in_between_2 = np.logical_and(in_between, np.logical_and(abias1 >= 1500, abias2 >= 1500))
    close_to_level_2 = np.logical_and(abias1 > 2000, abias2 < 1500)
    close_to_level_1 = abias1 < 1500
    lower = np.logical_and(abias2 > 1500, y[use_k] < x2[use_k])
    above = np.logical_and(abias1 > 1500, y[use_k] > x[use_k])
    print("between", np.sum(in_between)*100.0/len(in_between), np.sum(in_between_2)*100.0/len(in_between_2))
    print("det layer two", np.sum(close_to_level_2)*100.0/len(close_to_level_2))
    print("det layer one", np.sum(close_to_level_1)*100.0/len(close_to_level_1))
    print("lower", np.sum(lower)*100.0/len(lower))
    print("above", np.sum(above)*100.0/len(above))

    sort_ind = np.argsort(np.where(abias1 < abias2, abias1, abias2))
    print( np.mean(np.where(abias1<abias2, abias1, abias2)), np.mean(abias1), np.mean(abias2))
    print(  np.sum(np.where(abias1<abias2, abias1, abias2)<1000)*100.0/len(abias1))
    print(  np.sum(np.where(abias1<abias2, abias1, abias2)<1500)*100.0/len(abias1)) #(70%)

    sort_ind_top = np.argsort(abias1[closer_to_top])
    sort_ind_2 = np.argsort(abias2[closer_to_2])
    plt.plot(x[use_k][closer_to_top][sort_ind_top], 'g.')
    plt.plot(x2[use_k][closer_to_top][sort_ind_top], 'r.')
    plt.plot(y[use_k][closer_to_top][sort_ind_top], 'b.')
    ax = fig.add_subplot(322)
    plt.plot(x[use_k][closer_to_2][sort_ind_2], 'g.')
    plt.plot(x2[use_k][closer_to_2][sort_ind_2], 'r.')
    plt.plot(y[use_k][closer_to_2][sort_ind_2], 'b.')
    ax = fig.add_subplot(323)

    plt.plot(dist[closer_to_top][sort_ind_top], '.c')
    plt.plot(abias1[closer_to_top][sort_ind_top], 'k')
    ax = fig.add_subplot(324)

    plt.plot(dist[closer_to_2][sort_ind_2], '.c')
    plt.plot(abias2[closer_to_2][sort_ind_2], 'k')
    # plt.plot(abias1[sort_ind], 'g.')
    # plt.plot(abias2[sort_ind], 'r.')
    # plt.plot(np.where(abias1<abias2, abias1, abias2)[sort_ind], 'k')
    # plt.show()


def my_adjust_axis(ax):
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.grid(True)
    leg = plt.legend(loc="upper right", markerscale=1., numpoints=1, scatterpoints=1, bbox_to_anchor=(1.2, 1.05), framealpha=1.0, frameon=True)
    leg.get_frame().set_edgecolor('w')
    leg.get_frame().set_facecolor('w')
    # leg.get_frame().set_linewidth(0.0)
    ax.set_ylim(0, 10)
    ax.set_xlim(-4, 6)
    plt.yticks(np.arange(0, 9, 2.0))
    plt.yticks(np.arange(0, 9, 2.0))


def my_legend_text_6(data, text = "text"):
    # label1 = text
    label2 = "bias={:d}m\n".format(np.int(np.mean(data)))
    label3 = "STD={:d}m\n".format(np.int(np.std(data)))
    label4 = "MAE={:d}m\nIQR={:d}m\nQ2={:d}m\nPE0.5={:d}".format(
        # np.int(my_rms(data)),
        np.int(my_mae(data)),
        np.int(my_iqr(data)),
        np.int(np.median(data)),
        np.int(my_pe500m(data))
    )
    label = label2 + label3 + label4 +'%'
    return label


def my_make_plot_example(bias, use, label_str):
    plt.style.use('seaborn-white')
    fig = plt.figure(figsize=(11, 9))
    # plt.suptitle("CTTH error distributions not well described by RMS and bias")
    n_pix = 1000000
    ax = fig.add_subplot(221)
    temp_data = np.random.normal(900, 1600, n_pix)
    temp_data2 = np.concatenate([np.random.normal(-900, 200, int(0.5*n_pix)), np.random.normal( + 900, 200, int(0.5*n_pix))])
    hist_heights, x_m, dummy = my_hist(temp_data, None, bmin= -20*1000, bmax=20*1000, delta_h=100.0)
    hist_heights2, x_m, dummy = my_hist(temp_data2, None, bmin= -20*1000, bmax=20*1000, delta_h=100.0)
    hist_heights3, x_m, dummy = my_hist(bias, use, bmin= -20*1000, bmax=20*1000, delta_h=100.0)
    x_ = x_m*0.001
    ax = fig.add_subplot(221)
    # S1
    ax.set_title('a) Within threshold accuracy', loc='left')
    ax.fill(x_, hist_heights, color='silver',
             label = my_legend_text_6(temp_data, "Gaussian"))


#   plt.plot([0.001*np.mean(temp_data), 0.001*np.mean(temp_data)], [0, 2.4], 'k:')
    ax.set_ylabel('Percent')
    my_adjust_axis(ax)
    ax.set_xticklabels([])
    # S2
    ax = fig.add_subplot(222)
    ax.set_title('b) Within target accuracy', loc='left')
    ax.fill(x_, hist_heights2, color='grey',
            label = my_legend_text_6(temp_data2, "Bi-modal"))
    my_adjust_axis(ax)
    ax.set_yticklabels([])
    ax.set_xticklabels([])
    # S3
    ax = fig.add_subplot(223)
    ax.set_title('c) Outside threshold accuracy', loc='left')
    bias_i = bias[use]
    plt.plot(x_, hist_heights3, "r-",
             label = my_legend_text_6(bias_i, "PPS S-NPP"))
    # plt.plot([0.001*np.mean(bias_i), 0.001*np.mean(bias_i)], [0, 2], 'r:')
    my_adjust_axis(ax)
    ax.set_ylabel('Percent')
    ax.set_xlabel('error (km)')
    # S4
    ax = fig.add_subplot(224)
    ax.set_title('d) The worst (in STD) is the best', loc='left')
    ax.fill(x_, hist_heights, color='silver', label='Gaussian')
    ax.fill(x_, hist_heights2, color='grey', label='Bimodal')
    plt.plot(x_, hist_heights3, "r-", label='NN-CTTH')
    ax.set_xlabel('error (km)')
    my_adjust_axis(ax)
    ax.set_yticklabels([])
    plt.savefig(ADIR + "/PICTURES_FROM_PYTHON/VAL_2018_PLOTS/val_report_ctth_error_dist_%s.png"%(label_str), bbox_inches='tight')
    plt.savefig(ADIR + "/PICTURES_FROM_PYTHON/VAL_2018_PLOTS/val_report_ctth_error_dist_%s.pdf"%(label_str), bbox_inches='tight')
    # plt.show()


def my_make_plot_example_aprox(bias, use, label_str, match_calipso):
    def my_legend_text(data, text = "text"):
        # label1 = text
        label2 = "bias: {:d}\n".format(np.int(np.mean(data)))
        label3 = "STD: {:d}\n".format(np.int(np.std(data)))
        label4 = "IQR: {:d}\nQ2: {:d}".format(
            # np.int(my_rms(data)),
            # np.int(my_mae(data)),
            np.int(my_iqr(data)),
            np.int(np.median(data)),
            # np.int(my_pe500m(data))
        )
        label = label2 + label3 + label4
        return label
    bias_i = bias[use]
    n_pix = 1000000
    N_all = np.sum(use)
    temp_data_gs = []
    temp_data_iqrs = []
    from utils.get_flag_info import get_calipso_clouds_of_type_i
    for type_i in range(0, 8):
        # if type_i ==1:
        #   continue
        is_type_i = np.logical_and(use, get_calipso_clouds_of_type_i(match_calipso, calipso_cloudtype=type_i))
        n_pix_i = np.int(np.sum(is_type_i)*1.0/N_all*n_pix)
        temp_data_gi = np.random.normal(np.mean(bias[is_type_i]), np.std(bias[is_type_i]), n_pix_i)
        temp_data_iqri = np.random.normal(np.median(bias[is_type_i]), my_iqr(bias[is_type_i])*20.0/27, n_pix_i)
        if len(temp_data_gs) == 0:
            temp_data_gs = temp_data_gi
            temp_data_iqrs = temp_data_iqri
        else:
            temp_data_gs = np.concatenate([temp_data_gi, temp_data_gs ])
            temp_data_iqrs = np.concatenate([temp_data_iqri, temp_data_iqrs])

    temp_data_g = np.random.normal(np.mean(bias_i), np.std(bias_i), n_pix)
    temp_data_iqr = np.random.normal(np.median(bias_i), my_iqr(bias_i)*20.0/27, n_pix)

    # temp_data2 = np.concatenate([np.random.normal(-900, 200, int(0.5*n_pix)), np.random.normal(+900, 200, int(0.5*n_pix))])
    hist_heights_pps, x_m, dummy = my_hist(bias_i, None, bmin= -20*1000, bmax=20*1000, delta_h=100.0)

    hist_heights_g, x_m, dummy = my_hist(temp_data_g, None, bmin= -20*1000, bmax=20*1000, delta_h=100.0)
    hist_heights_iqr, x_m, dummy = my_hist(temp_data_iqr, None, bmin= -20*1000, bmax=20*1000, delta_h=100.0)
    hist_heights_gs, x_m, dummy = my_hist(temp_data_gs, None, bmin= -20*1000, bmax=20*1000, delta_h=100.0)
    hist_heights_iqrs, x_m, dummy = my_hist(temp_data_iqrs, None, bmin= -20*1000, bmax=20*1000, delta_h=100.0)
    # hist_heights3, x_m, dummy = my_hist(bias, use, bmin=-20*1000, bmax=20*1000, delta_h=100.0)
    x_ = x_m*0.001

    def gaussian_fit(x, y, bias):
        n = len(x)                          #the number of data
        mean = np.median(bias)                   #note this correction
        sigma = my_iqr(bias)*21.0/27        #note this correction
        from scipy.optimize import curve_fit
        from scipy import asarray as ar, exp

        def gaus(x, a, x0, sigma):
            return 8*exp( - np.abs(x - x0)*np.sqrt(2)/(sigma))

        # popt, pcov = curve_fit(gaus, x, y, p0=[1, mean, sigma])

        plt.plot(x, y, 'b+:', label='data')
        plt.plot(x, gaus(x, 1, mean, sigma), 'ro:', label='fit')
        plt.show()
        return 1, mean, sigma

    dummy, mean_g, sigma_g = gaussian_fit(x_m, hist_heights_pps, bias_i)
    temp_data_fitted = np.random.laplace(np.median(bias), my_iqr(bias)*32/26.0, n_pix)
    hist_heights_fitted, x_m, dummy = my_hist(temp_data_fitted, None, bmin= -20*1000, bmax=20*1000, delta_h=100.0)

    fig = plt.figure(figsize=(7, 5))
    ax = fig.add_subplot(111)
    ax.fill(x_, hist_heights_fitted, color='silver',
            label = "\nGaussian\n" + my_legend_text_6(temp_data_fitted, ))
    plt.plot(x_, hist_heights_pps, "r-",
             label = "NN-CTTH\n"+my_legend_text_6(bias_i, "PPS"  ))
    ax.set_xlim(-4, 6)
    leg = plt.legend(loc="upper right", markerscale=2., numpoints=1, scatterpoints=1, bbox_to_anchor=(1.13, 1.2), framealpha=1.0, frameon=True)
    leg.get_frame().set_facecolor('w')
    plt.savefig(ADIR + "/PICTURES_FROM_PYTHON/VAL_2018_PLOTS/val_report_ctth_error_and_gaussian_fitted_%s.png"%(label_str), bbox_inches='tight')
    plt.savefig(ADIR + "/PICTURES_FROM_PYTHON/VAL_2018_PLOTS/val_report_ctth_error_and_gaussian_fitted_%s.pdf"%(label_str), bbox_inches='tight')

    plt.style.use('seaborn-white')
    fig = plt.figure(figsize=(7, 5))
    # plt.suptitle("CTTH error distributions not well described by RMS and bias")
    ax = fig.add_subplot(111)
    # S1
    # ax.set_title('PPS-v2018', loc='left', zorder=30)
    ax.fill(x_, hist_heights_g, color='silver',
             label = "\nGaussian\n" + my_legend_text(temp_data_g, ))
    plt.plot(x_, hist_heights_pps, "r-",
             label = "NN-CTTH\n"+my_legend_text(bias_i, "PPS"  ))
    plt.plot([0.001*np.mean(bias_i), 0.001*np.mean(bias_i)], [0, 2], 'r:')
    ax.set_xlim(-4, 6)
    ax.set_ylabel('Percent')
    ax.set_xlabel('Error (km)')
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.grid(True)
    leg = plt.legend(loc="upper right", markerscale=2., numpoints=1, scatterpoints=1, bbox_to_anchor=(1.13, 1.2), framealpha=1.0, frameon=True)
    leg.get_frame().set_facecolor('w')
    leg.get_frame().set_linewidth(0.0)
    # import pdb
    # pdb.set_trace()
    print("inside {:3.1f}".format(np.sum(hist_heights_pps[np.logical_and(x_>=-4, x_<=6)])))
    print("gaussian inside {:3.1f}".format(np.sum(hist_heights_g[np.logical_and(x_>=-4, x_<=6)])))
    # plt.legend(frameon=False)

    ax.set_ylim(0, 10)
    ax.set_xlim(-4, 6)
    plt.yticks(np.arange(0, 9, 2.0))
    plt.savefig(ADIR + "/PICTURES_FROM_PYTHON/VAL_2018_PLOTS/val_report_ctth_error_and_gaussian_%s.png"%(label_str), bbox_inches='tight')
    plt.savefig(ADIR + "/PICTURES_FROM_PYTHON/VAL_2018_PLOTS/val_report_ctth_error_and_gaussian_%s.pdf"%(label_str), bbox_inches='tight')

    plt.style.use('seaborn-white')
    fig = plt.figure(figsize=(12, 4.5))
    # plt.suptitle("CTTH error distributions not well described by RMS and bias")
    ax = fig.add_subplot(121)
    # S1
    ax.grid(True)
    ax.set_title('PPS-v2018', loc='left', zorder=30)
    # ax.fill(x_, hist_heights_g, color='silver',
    #        label = my_legend_text(temp_data_g, "Gaussian (STD)"))
    plt.fill(x_, hist_heights_pps, "r-",  #alpha = 0.5,
             label = my_legend_text_6(bias_i, "PPS"  ))
    plt.plot(x_-0.001*np.mean(bias_i), hist_heights_pps, "b-")
    ax.set_xlim(-4, 6)
    ax.set_ylabel('Percent')
    ax.set_xlabel('Error (km)')
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)

    leg = plt.legend(loc="upper right", markerscale=2., numpoints=1, scatterpoints=1, bbox_to_anchor=(1.13, 1.2), framealpha=1.0, frameon=True)
    leg.get_frame().set_facecolor('w')
    leg.get_frame().set_linewidth(0.0)
    ax.set_ylim(0, 10)
    ax.set_xlim(-4, 6)
    plt.yticks(np.arange(0, 9, 2.0))
    ax = fig.add_subplot(122)
    ax.grid(True)
    # S1
    ax.set_title('Bias "corrected"', loc='left', zorder=30)
    # ax.fill(x_-0.001*np.mean(temp_data_g), hist_heights_g, color='silver',
    #        label = my_legend_text(temp_data_g-0.001*np.mean(temp_data_g), "Gaussian (STD)"))
    plt.fill(x_-0.001*np.mean(bias_i), hist_heights_pps, "b-", #alpha = 0.5,
             label = my_legend_text_6(bias_i-np.mean(bias_i), "PPS"  ))
    plt.plot(x_, hist_heights_pps, "r-")
    ax.set_xlim(-4, 6)
    # ax.set_ylabel('Percent')
    ax.set_xlabel('Error (km)')
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    leg = plt.legend(loc="upper right", markerscale=2., numpoints=1, scatterpoints=1, bbox_to_anchor=(1.13, 1.2), framealpha=1.0, frameon=True)
    leg.get_frame().set_facecolor('w')
    leg.get_frame().set_linewidth(0.0)
    ax.set_ylim(0, 10)
    ax.set_xlim(-4, 6)
    plt.yticks(np.arange(0, 9, 2.0))
    plt.savefig(ADIR + "/PICTURES_FROM_PYTHON/VAL_2018_PLOTS/val_report_ctth_error_and_gaussian_and_biascorr_%s.png"%(label_str), bbox_inches='tight')
    plt.savefig(ADIR + "/PICTURES_FROM_PYTHON/VAL_2018_PLOTS/val_report_ctth_error_and_gaussian_and_biascorr_%s.pdf"%(label_str), bbox_inches='tight')

    plt.close()
    plt.style.use('seaborn-white')
    fig = plt.figure(figsize=(15, 11))
    # plt.suptitle("CTTH error distributions not well described by RMS and bias")
    ax = fig.add_subplot(221)
    # S1
    ax.set_title('a) Equal bias/std')
    ax.fill(x_, hist_heights_g, color='silver',
             label = my_legend_text(temp_data_g, "Gaussian (STD)"))
    plt.plot(x_, hist_heights_pps, "r-",
             label = "\n" + my_legend_text(bias_i, "PPS"))
    ax.set_ylabel('Percent')
    my_adjust_axis(ax)
    ax.set_xticklabels([])
    # S2
    ax = fig.add_subplot(222)
    ax.set_title('b) Equal IQR/median')
    ax.fill(x_, hist_heights_iqr, color='grey',
            label = "\n" +my_legend_text(temp_data_iqr, "Gaussian (IQR)"))
    plt.plot(x_, hist_heights_pps, "r-")#, label='PPS')
    my_adjust_axis(ax)
    ax.set_yticklabels([])
    ax.set_xticklabels([])
    # S3
    ax = fig.add_subplot(223)
    ax.set_title('c) Equal bias/std, sum of types')

    plt.fill(x_, hist_heights_gs, color='silver',
             label = "\n" + my_legend_text(temp_data_gs, "Gaussian \Sum (STD)"))
    plt.plot(x_, hist_heights_pps, "r-")#, label='PPS')
    my_adjust_axis(ax)
    ax.set_ylabel('Percent')
    ax.set_xlabel('error (km)')
    # S4
    ax = fig.add_subplot(224)
    ax.set_title('d) Equal IQR/median, sum over cloud types')
    ax.fill(x_, hist_heights_iqrs, color='grey',
            label = "\n" + my_legend_text(temp_data_iqrs, "Gaussian \Sum (IQR)"))
    plt.plot(x_, hist_heights_pps, "r-")#, label='PPS')

    ax.set_xlabel('error (km)')
    my_adjust_axis(ax)
    ax.set_yticklabels([])
    plt.savefig(ADIR + "/PICTURES_FROM_PYTHON/VAL_2018_PLOTS/val_report_ctth_error_aprox_dist_%s.png"%(label_str), bbox_inches='tight')
    plt.savefig(ADIR + "/PICTURES_FROM_PYTHON/VAL_2018_PLOTS/val_report_ctth_error_aprox_dist_%s.pdf"%(label_str), bbox_inches='tight')
    # plt.show()


def my_print_one_line(out_file_h, bias, x, y, use_ind, compare_name, flag_key):
    out_line = "%s_%s %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %3.1f "%(
        compare_name, flag_key,
        np.mean(bias[use_ind]),
        np.median(bias[use_ind]),
        my_iqr(bias[use_ind]),
        my_pe500m(bias[use_ind]),
        my_pe1000m(bias[use_ind]),
        my_mae(bias[use_ind]),
        np.std(bias[use_ind]),
        len(bias[use_ind]),
        my_rms(bias[use_ind]),
        np.mean(x[use_ind]),
        np.mean(y[use_ind]),
        my_pe250m(bias[use_ind]),
        my_pe2000m(bias[use_ind]),
        my_pe2500m(bias[use_ind]),
        my_pe5000m(bias[use_ind]),
        half_sample_mode(bias[use_ind]),
        skew(bias[use_ind]), #,
        # kurtosis(bias[use_ind])
    )
    out_file_h.write(out_line)
    out_file_h.write("\n")
    print(out_line)


def print_all_cloudsat(match_obj, compare, compare_name = "unknown"):
    from utils.get_flag_info import get_cloudsat_low_medium_high_classification
    x = match_obj.cloudsat.all_arrays['validation_height']
    y = match_obj.imager.all_arrays['imager_ctth_m_above_seasurface']
    mhl = get_cloudsat_low_medium_high_classification(match_obj)
    use = np.logical_and(x >= 0, np.logical_and(y > - 9, y < 65000))
    if mhl is None:
        mhl = {}
    mhl["all"] = use.copy()
    bias = y - x
    out_file_h.write(compare_name + " CPR (CloudSat) :\n")
    for flag_key in ["all", "low_clouds", "medium_clouds", "high_clouds"]:
        if flag_key not in mhl.keys():
            continue
        else:
            use_i = np.logical_and(use, np.logical_and(mhl[flag_key], use))
            my_print_one_line(out_file_h, bias, x, y, use_i, "all", flag_key)


def print_all(match_obj, compare, compare_name = "unknown"):
    # x = getattr(plt_obj, truth)
    # y = getattr(plt_obj, compare)

    x = match_obj.calipso.all_arrays['validation_height']
    x2 = match_obj.calipso.all_arrays['layer_top_altitude'][:, 1]*1000
    y = match_obj.imager.all_arrays['imager_ctth_m_above_seasurface']
    print( np.max(y), np.max(y[y<65000]))
    # pressure_c = match_obj.calipso.all_arrays['layer_top_pressure'][:, 0]
    low_clouds = get_calipso_low_clouds(match_obj)
    high_clouds = get_calipso_high_clouds(match_obj)
    medium_clouds = get_calipso_medium_clouds(match_obj)
    mhl = get_calipso_low_medium_high_classification(match_obj)

    use = np.logical_and(x >= 0, np.logical_and(y > - 9, y < 65000))
    use = np.logical_and(use, np.not_equal(match_obj.calipso.all_arrays['feature_classification_flags'][:, 0], 1))
    mhl["all"] =use
    use_inversion = get_inversion_info_pps2014(match_obj.imager.all_arrays["cloudtype_status"])

    mhl["high_clouds_tp_thin"] = np.logical_and(
        mhl["high_clouds_tp"],
        np.logical_and(match_obj.calipso.all_arrays['feature_optical_depth_532_top_layer_5km']>=0,
                       match_obj.calipso.all_arrays['feature_optical_depth_532_top_layer_5km'] <= 0.225))
    mhl["high_clouds_tp_not_thin"] = np.logical_and(
        mhl["high_clouds_tp"], ~mhl["high_clouds_tp_thin"])

    mhl["medium_clouds_tp_thin"] = np.logical_and(
        mhl["medium_clouds_tp"],
        np.logical_and(match_obj.calipso.all_arrays['feature_optical_depth_532_top_layer_5km']>=0,
                      match_obj.calipso.all_arrays['feature_optical_depth_532_top_layer_5km'] <= 0.225))
    mhl["medium_clouds_tp_not_thin"] = np.logical_and(
        mhl["medium_clouds_tp"], ~mhl["medium_clouds_tp_thin"])

    mhl["low_clouds_tp_thin"] = np.logical_and(
        mhl["low_clouds_tp"],
        np.logical_and(match_obj.calipso.all_arrays['feature_optical_depth_532_top_layer_5km']>=0,
                      match_obj.calipso.all_arrays['feature_optical_depth_532_top_layer_5km'] <= 0.225))
    mhl["low_clouds_tp_not_thin"] = np.logical_and(
        mhl["low_clouds_tp"], ~mhl["low_clouds_tp_thin"] )

    mhl["all_clouds_tp_thin"] = np.logical_and(
        mhl["clouds_tp"],
        np.logical_and(match_obj.calipso.all_arrays['feature_optical_depth_532_top_layer_5km']>=0,
                       match_obj.calipso.all_arrays['feature_optical_depth_532_top_layer_5km'] <= 0.225))
    mhl["all_clouds_tp_not_thin"] = np.logical_and(
        mhl["clouds_tp"], ~mhl["all_clouds_tp_thin"] )

    mhl["all_clouds_tp_bad2"] = np.logical_and(
        mhl["high_clouds_tp_thin"],
        match_obj.calipso.all_arrays['number_layers_found']>=2)

    mhl["all_clouds_tp_bad3"] = np.logical_and(
        mhl["high_clouds_tp_thin"],
        match_obj.calipso.all_arrays['number_layers_found']==1)

    use_low = np.logical_and(use, low_clouds)
    use_medium = np.logical_and(use, medium_clouds)
    use_high = np.logical_and(use, high_clouds)
    bias = y - x#+1465
    bias_second_layer = y - x2
    abias = np.abs(bias)

    # my_make_plot_example(bias, use, compare_name)
    # my_make_plot2(y, x, x2, mhl, mhl["all"])
    my_make_plot_example_aprox(bias, use, compare_name, match_obj)

    from scipy import  ndimage
    maxct = ndimage.filters.maximum_filter1d(match_obj.imager.cloudtype, size=9)
    minct = ndimage.filters.minimum_filter1d(match_obj.imager.cloudtype, size=9)
    val_geo = np.equal(maxct, minct)
    # match_obj.calipso.layer_top_pressure[:, 0][match_obj.calipso.layer_top_pressure[:, 0]<0] =1200
    # match_obj.calipso.layer_top_altitude[:, 0][match_obj.calipso.layer_top_altitude[:, 0]<0] =0
    if hasattr(match_obj, 'calipso'):
        var_pressure = (ndimage.filters.maximum_filter1d(match_obj.calipso.layer_top_pressure[:, 0], size=9) -
                        ndimage.filters.minimum_filter1d(match_obj.calipso.layer_top_pressure[:, 0], size=9))
        val_geo = np.logical_and(
            val_geo,
            var_pressure < 200) #Pressure variation less than 200hPa
    var_pressure = (ndimage.filters.maximum_filter1d(match_obj.calipso.layer_top_pressure[:, 0], size=9) -
                    ndimage.filters.minimum_filter1d(match_obj.calipso.layer_top_pressure[:, 0], size=9))
    var_height = (ndimage.filters.maximum_filter1d(match_obj.calipso.layer_top_altitude[:, 0]*1000, size=9) -
                  ndimage.filters.minimum_filter1d(match_obj.calipso.layer_top_altitude[:, 0]*1000, size=9))
    val_geo2 = var_pressure < 100
    sunz = np.array(match_obj.imager.all_arrays['sunz'])

    """
    fig = plt.figure(figsize=(15, 11))
    print match_obj.calipso.all_arrays['feature_classification_flags'][use][:, 0]
    cflag_full = match_obj.calipso.all_arrays['feature_classification_flags'][:, 0]
    cflag = match_obj.calipso.all_arrays['feature_classification_flags'][:, 0][use]
    feature_array = (4*np.bitwise_and(np.right_shift(cflag, 11), 1) +
                     2*np.bitwise_and(np.right_shift(cflag, 10), 1) +
                     1*np.bitwise_and(np.right_shift(cflag, 9), 1))
    feature_array = (4*np.bitwise_and(np.right_shift(cflag, 15), 1) +
                     2*np.bitwise_and(np.right_shift(cflag, 14), 1) +
                     1*np.bitwise_and(np.right_shift(cflag, 13), 1))
    feature_array_full = (2*np.bitwise_and(np.right_shift(cflag_full, 4), 1) +
                          1*np.bitwise_and(np.right_shift(cflag_full, 3), 1))
    feature_array = (2*np.bitwise_and(np.right_shift(cflag, 4), 1) +
                     1*np.bitwise_and(np.right_shift(cflag, 3), 1))
    print np.mean(abias[use][feature_array==0]), len(abias[use][feature_array==0])
    print np.mean(abias[use][feature_array==1]), len(abias[use][feature_array==1])
    print np.mean(abias[use][feature_array==2]), len(abias[use][feature_array==2])
    print np.mean(abias[use][feature_array==3]), len(abias[use][feature_array==3])
    print my_rms(abias[use][feature_array==0]), len(abias[use][feature_array==0])
    print my_rms(abias[use][feature_array==1]), len(abias[use][feature_array==1])
    print my_rms(abias[use][feature_array==2]), len(abias[use][feature_array==2])
    print my_rms(abias[use][feature_array==3]), len(abias[use][feature_array==3])

    fig = plt.figure(figsize=(15, 11))
    ax = fig.add_subplot(111)
    # plt.plot(var_pressure[use], abias[use], 'b.', alpha=0.02)
    # plt.plot(var_height[use], abias[use], 'b.', alpha=0.02)
    # plt.plot(match_obj.calipso.all_arrays['feature_optical_depth_532_top_layer_5km'][use], abias[use], 'b.', alpha=0.02)
    # plt.plot(match_obj.calipso.all_arrays['cfc_single_shots_1km_from_5km_file'][use]+0.01), abias[use], 'b.', alpha=0.02)
    # ax.set_xlim([-1.0, 1.0])
    plt.plot(feature_array, abias[use], 'b.', alpha=0.02)
    # plt.plot(match_obj.calipso.layer_top_altitude[:, 0][use], abias[use], 'r.', alpha=0.05)
    plt.show()
    """
    # logger.info("Getting day/night info from sunz")
    if np.max(sunz) < 20:
        sunz =sunz*100.0
    day_flag = np.where(np.less_equal(sunz, 80), 1, 0)
    night_flag = np.where(np.greater_equal(sunz, 95), 1, 0)
    twilight_flag = np.where(
        np.logical_and(np.greater(sunz, 80),
                       np.less(sunz, 95)),
        1, 0)
    out_file_h.write(compare_name + " CALIOP:\n")
    # use = np.logical_and(use, feature_array_full==3)
    my_list = sorted(mhl.keys())
    my_list = ["clouds_tp", "low_clouds_tp", "medium_clouds_tp", "high_clouds_tp",
               # "all_clouds_tp_not_thin", "low_clouds_tp_not_thin", "medium_clouds_tp_not_thin", "high_clouds_tp_not_thin",
               # "all_clouds_tp_thin", "low_clouds_tp_thin", "medium_clouds_tp_thin", "high_clouds_tp_thin",
               "clouds_op", "low_clouds_op", "medium_clouds_op", "high_clouds_op",
               "all_clouds_tp_bad2", "all_clouds_tp_bad3"]
    for flag, compare_name in zip( [use, day_flag, night_flag, twilight_flag, val_geo, val_geo2],
                                   ["all", "day", "night", "twilight", "geo-style", "no-edges"]):
        use_i = np.logical_and(use, flag)
        my_print_one_line(out_file_h, bias, x, y, use_i, compare_name, "")
    for flag_key in my_list:
        use_i = np.logical_and(use, np.logical_and(mhl[flag_key], use))
        my_print_one_line(out_file_h, bias, x, y, use_i, "all", flag_key)
    for flag_key in my_list:
        use_i = np.logical_and(use, np.logical_and(mhl[flag_key], use_inversion))
        my_print_one_line(out_file_h, bias, x, y, use_i, "all", flag_key+'inversion')
    use_i = np.logical_and(use, np.logical_and(mhl["all_clouds_tp_bad2"], use))
    my_print_one_line(out_file_h, bias_second_layer, x2, y, use_i, "all", "all_clouds_tp_bad2")
    plt.close('all')

if __name__ == "__main__":

    from matchobject_io import read_files
    BASE_DIR = ADIR + "/DATA_MISC/reshaped_files_validation_2018/"
    ROOT_DIR_v2014 = (BASE_DIR + "global_modis_v2014_created20180920/Reshaped_Files_merged_caliop/eos2/1km/2010/*/*%s*h5")
    ROOT_DIR_v2014_clsat = (BASE_DIR + "global_modis_v2014_created20180920/Reshaped_Files_merged_cloudsat/eos2/1km/2010/*/*%s*h5")
    ROOT_DIR_v2018 = (BASE_DIR + "global_modis_v2018_created20180920/Reshaped_Files_merged_caliop/eos2/1km/2010/*/*%s*h5")
    ROOT_DIR_v2018_clsat = (BASE_DIR + "global_modis_v2018_created20180920/Reshaped_Files_merged_cloudsat/eos2/1km/2010/*/*%s*h5")

    ROOT_DIR_v2014_NPP = (BASE_DIR + "global_viirs_v2014_created20180914/Reshaped_Files_merged_caliop/npp/1km/2015/*/*h5")
    ROOT_DIR_v2018_NPP = (BASE_DIR + "global_viirs_v2018_created20180907/Reshaped_Files_merged_caliop/npp/1km/2015/*/*h5")
    ROOT_DIR_v2014_NPP_clsat = (BASE_DIR + "global_viirs_v2014_created20180914/Reshaped_Files_merged_cloudsat/npp/1km/2015/*/*h5")
    ROOT_DIR_v2018_NPP_clsat = (BASE_DIR + "global_viirs_v2018_created20181002_new_cmaprobv5/Reshaped_Files_merged_cloudsat/npp/1km/2015/*/*h5")
    ROOT_DIR_v2014_GAC = (BASE_DIR + "global_gac_v2014_created20180927/Reshaped_Files/noaa18/5km/%s/*cali*h5")
    ROOT_DIR_v2018_GAC = (BASE_DIR + "global_gac_v2018_created20180927/Reshaped_Files/noaa18/5km/%s/*cali*h5")
    ROOT_DIR_v2014_GAC_clsat = (BASE_DIR + "global_gac_v2014_created20180927/Reshaped_Files/noaa18/5km/2009/*clouds*h5")
    ROOT_DIR_v2018_GAC_clsat = (BASE_DIR + "global_gac_v2018_created20180927/Reshaped_Files/noaa18/5km/2009/*clouds*h5")
    ROOT_DIR = ROOT_DIR_v2018_clsat

# between 27.16770203522559 12.67069745562212
# between 38.70328344488811 22.066300881905146

    files = glob(ROOT_DIR_v2014_NPP)
    out_file_h.write("NPP-v2014\n")
    match_obj = read_files(files)
    print_all(match_obj, None, "NPPv2014")

    files = glob(ROOT_DIR_v2018_NPP)
    out_file_h.write("NPP-v2018\n")
    match_obj = read_files(files)
    print_all(match_obj, None, "NPPv2018")
    b=a

    files = glob(ROOT_DIR%("20100201"))
    files = files + glob(ROOT_DIR%("20100401"))
    files = files + glob(ROOT_DIR%("20100601"))
    files = files + glob(ROOT_DIR%("20100801"))
    files = files + glob(ROOT_DIR%("20101001"))
    files = files + glob(ROOT_DIR%("20101201"))
    out_file_h.write("MODIS-C6\n")
    match_obj = read_files(files, truth='cloudsat')
    match_obj.imager.all_arrays['imager_ctth_m_above_seasurface'] = match_obj.modis.all_arrays["height"]
    print_all_cloudsat(match_obj, None, "MODIS-C6")

    ROOT_DIR = ROOT_DIR_v2018
    files = glob(ROOT_DIR%("20100201"))
    files = files + glob(ROOT_DIR%("20100401"))
    files = files + glob(ROOT_DIR%("20100601"))
    files = files + glob(ROOT_DIR%("20100801"))
    files = files + glob(ROOT_DIR%("20101001"))
    files = files + glob(ROOT_DIR%("20101201"))
    out_file_h.write("MODIS-C6\n")
    match_obj = read_files(files)
    match_obj.imager.all_arrays['imager_ctth_m_above_seasurface'] = match_obj.modis.all_arrays["height"]
    print_all(match_obj, None, "MODIS-C6")

    files = glob(ROOT_DIR_v2014_GAC%("2006"))
    files = files + glob(ROOT_DIR_v2014_GAC%("2009"))
    out_file_h.write("GAc-v2014\n")
    match_obj = read_files(files)
    print_all(match_obj, None, "GACv2014")
    files = glob(ROOT_DIR_v2018_GAC%("2006"))
    files = files + glob(ROOT_DIR_v2018_GAC%("2009"))
    out_file_h.write("GAC-v2018\n")
    match_obj = read_files(files)
    print_all(match_obj, None, "GACv2018")

    files = glob(ROOT_DIR_v2014_GAC_clsat) #only 2009
    out_file_h.write("GAc-v2014\n")
    match_obj = read_files(files, truth='cloudsat')
    print_all_cloudsat(match_obj, None, "GACv2014")
    files = glob(ROOT_DIR_v2018_GAC_clsat) #only 2009
    out_file_h.write("GAC-v2018\n")
    match_obj = read_files(files, truth='cloudsat')
    print_all_cloudsat(match_obj, None, "GACv2018")

    files = glob(ROOT_DIR_v2018_NPP_clsat)
    out_file_h.write("NPP-v2018\n")
    match_obj = read_files(files, truth='cloudsat')
    print_all_cloudsat(match_obj, None, "NPPv2018")
    files = glob(ROOT_DIR_v2014_NPP_clsat)
    out_file_h.write("NPP-v2014\n")
    match_obj = read_files(files, truth='cloudsat')
    print_all_cloudsat(match_obj, None, "NPPv2014")

    ROOT_DIR = ROOT_DIR_v2014_clsat
    files = glob(ROOT_DIR%("20100201"))
    files = files + glob(ROOT_DIR%("20100401"))
    files = files + glob(ROOT_DIR%("20100601"))
    files = files + glob(ROOT_DIR%("20100801"))
    files = files + glob(ROOT_DIR%("20101001"))
    files = files + glob(ROOT_DIR%("20101201"))
    out_file_h.write("MODIS-v2014\n")
    match_obj = read_files(files, truth='cloudsat')
    print_all_cloudsat(match_obj, None, "MODISv2014")

    ROOT_DIR = ROOT_DIR_v2018_clsat
    files = glob(ROOT_DIR%("20100201"))
    files = files + glob(ROOT_DIR%("20100401"))
    files = files + glob(ROOT_DIR%("20100601"))
    files = files + glob(ROOT_DIR%("20100801"))
    files = files + glob(ROOT_DIR%("20101001"))
    files = files + glob(ROOT_DIR%("20101201"))
    out_file_h.write("MODIS-v2018\n")
    match_obj = read_files(files, truth='cloudsat')
    print_all_cloudsat(match_obj, None, "MODISv2018")

    ROOT_DIR = ROOT_DIR_v2014
    files = glob(ROOT_DIR%("20100201"))
    files = files + glob(ROOT_DIR%("20100401"))
    files = files + glob(ROOT_DIR%("20100601"))
    files = files + glob(ROOT_DIR%("20100801"))
    files = files + glob(ROOT_DIR%("20101001"))
    files = files + glob(ROOT_DIR%("20101201"))
    out_file_h.write("MODIS-v2014\n")
    match_obj = read_files(files)
    print_all(match_obj, None, "eos2v2014")

    ROOT_DIR = ROOT_DIR_v2018
    files = glob(ROOT_DIR%("20100201"))
    files = files + glob(ROOT_DIR%("20100401"))
    files = files + glob(ROOT_DIR%("20100601"))
    files = files + glob(ROOT_DIR%("20100801"))
    files = files + glob(ROOT_DIR%("20101001"))
    files = files + glob(ROOT_DIR%("20101201"))
    out_file_h.write("MODIS-v2018\n")
    match_obj = read_files(files)
    print_all(match_obj, None, "eos2v2018")


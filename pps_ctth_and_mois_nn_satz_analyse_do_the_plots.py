import os
import re
import netCDF4
from glob import glob
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
matplotlib.rcParams.update(matplotlib.rcParamsDefault)
matplotlib.rcParams.update({'font.size': 18})
plt.rc('text', usetex=True)
plt.rc('font', family='serif')

from read_cloudproducts_and_nwp_pps import (read_ctth_nc, read_pps_angobj_nc)
import matplotlib
#matplotlib.use("TkAgg")

def read_satz_statistics_files(out_filenames):
    aggregated_data = {}
    for out_filename in out_filenames:
        print out_filename
        out_file_h = open(out_filename,'r')
        for line in out_file_h:
            #print line
            line = line.rstrip()
            line = line.lstrip()
            line_as_list = line.split(" ")
            #print   line_as_list
            start_with_number = re.match("[\d+\s+]+\d+$",line)
            data_line = re.match("^(\S+)\s+satz(\d+)_(\d\d)\s+([\d+\s+]+\d+)\s*$",line)
        
            if start_with_number:
                print " got the pressure array"
                pressure_plot  = np.array([np.float(i) for i in line_as_list])
                print pressure_plot 
            elif data_line:
                #print "got data", line
                ctth_type = data_line.group(1)
                if ctth_type not in aggregated_data.keys():
                    aggregated_data[ctth_type] = {}
                min_satz = np.int(data_line.group(2))
                max_satz = int(data_line.group(3))
                data =  np.array([np.float(i) for i in data_line.group(4).split(" ")])
                satz_step = max_satz-min_satz
                if min_satz not in aggregated_data[ctth_type].keys():
                    print min_satz
                    aggregated_data[ctth_type][min_satz] = data
                else:    
                    aggregated_data[ctth_type][min_satz] += data
    return aggregated_data, satz_step,pressure_plot          

def plot_pressure_local(out_filenames,PLOT_DIR):
    aggregated_data, satz_step,pressure_plot = read_satz_statistics_files(out_filenames)
    step = 50#hPa read this!!!
    #["CTTHold", "CTTHnnAvhrr", "CTTHnnAvhrrNoRTTOV"]):     
    fig = plt.figure(figsize = (12,12))  
    plt.suptitle('CTTH pressure dependence on satzenith angle')
    def plot_one_subplot(aggregated_data, satz_step,pressure_plot,  ctth_label, desc, label=False):
        ax.text(80,11, desc, bbox={'facecolor':'blue', 'alpha':0.1, 'pad':1})
        plot_00_20 = (aggregated_data[ctth_label][0] + 
                      aggregated_data[ctth_label][10])
        plot_20_40 = (aggregated_data[ctth_label][20] + 
                      aggregated_data[ctth_label][30])
        plot_40_60 = (aggregated_data[ctth_label][40] + 
                      aggregated_data[ctth_label][50])
        plot_60_80 = (aggregated_data[ctth_label][60] + 
                      aggregated_data[ctth_label][70])
        plot_00_20 = plot_00_20*100.0/np.sum(plot_00_20)
        plot_20_40 = plot_20_40*100.0/np.sum(plot_20_40)
        plot_40_60 = plot_40_60*100.0/np.sum(plot_40_60)
        plot_60_80 = plot_60_80*100.0/np.sum(plot_60_80)
        
        plt.plot(pressure_plot[:-1], plot_00_20, '-k.', label="satz: 0-20")
        plt.plot(pressure_plot[:-1], plot_20_40, '-ks', label="satz: 20-40")
        plt.plot(pressure_plot[:-1], plot_40_60, '-k*', label="satz: 40-60")
        plt.plot(pressure_plot[:-1], plot_60_80, '-ko', label="satz: 60-80")
        ax.set_ylabel('percent of ctth results')
        ax.set_ylim(0,15)
        ax.set_xlim(70,1200)
        if label:
            plt.legend(loc='upper right',bbox_to_anchor=(1.1, 1.05))
    ax = fig.add_subplot(221)
    plot_one_subplot(aggregated_data, satz_step, pressure_plot,"CTTHold", "CTTH old")
    ax = fig.add_subplot(222)
    plot_one_subplot(aggregated_data, satz_step, pressure_plot, "CTTHnnAvhrr", "CTTH nn", label=True)
    ax = fig.add_subplot(223)
    plot_one_subplot(aggregated_data, satz_step, pressure_plot,"CTTHnnAvhrrNoRTTOV", "CTTH nn no THR")
    ax = fig.add_subplot(224)
    plot_one_subplot(aggregated_data, satz_step, pressure_plot,"MODIS-C6", "MODIS-C6")

    ax.set_xlabel('pressure (hPa)')
    plt.savefig(PLOT_DIR + "ctth_satz_pressure_%s_second_try_.png"%('modis'), bbox_inches='tight')

def plot_pressure_modis_01(out_filenames,PLOT_DIR,month):
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')
    aggregated_data, satz_step,pressure_plot = read_satz_statistics_files(out_filenames)
    step = 50#hPa read this!!!
    #["CTTHold", "CTTHnnAvhrr", "CTTHnnAvhrrNoRTTOV"]):     
    fig = plt.figure(figsize = (19,10))  
    plt.suptitle('CTTH pressure dependence on satzenith angle', fontsize=22)
    def plot_one_subplot(aggregated_data, satz_step,pressure_plot,  ctth_label, desc, legend=False):
        ax.text(300,11.5, desc, bbox={'facecolor':'blue', 'alpha':0.0, 'pad':1})
        #plt.text(0.02, 0.90, desc, fontsize=12,transform=ax.transAxes, bbox=dict(facecolor='w', edgecolor='w', alpha=1.0))
        ax.grid(True)
        plot_a = (aggregated_data[ctth_label][0] ) 
                      #aggregated_data[ctth_label][10])
        plot_b = (aggregated_data[ctth_label][20] + 
                      aggregated_data[ctth_label][10])
        plot_c = (aggregated_data[ctth_label][40] + 
                      aggregated_data[ctth_label][30])
        plot_d = (aggregated_data[ctth_label][60] + 
                      aggregated_data[ctth_label][50])
        plot_calipso = aggregated_data["CALIPSO"][0]
        plot_a = plot_a*100.0/np.sum(plot_a)
        plot_b = plot_b*100.0/np.sum(plot_b)
        plot_c = plot_c*100.0/np.sum(plot_c)
        plot_d = plot_d*100.0/np.sum(plot_d)
        plot_calipso = plot_calipso*100.0/np.sum(plot_calipso) 
        plt.plot(pressure_plot[:-1], plot_calipso, 'b', alpha=0.2, label="CALIPSO") 
        plt.plot(pressure_plot[:-1], plot_a, '-k.', label="satz: 0-10")
        plt.plot(pressure_plot[:-1], plot_b, '-ks', label="satz: 10-30")
        plt.plot(pressure_plot[:-1], plot_c, '-k*', label="satz: 30-50")
        plt.plot(pressure_plot[:-1], plot_d, '-ko', label="satz: 50-70")
        #ax.set_ylabel('percent of ctth results')
        ax.set_ylim(0,13)
        ax.set_xlim(70,1200)
        if legend:
            plt.legend(loc='upper right',bbox_to_anchor=(1.3, 1.3),fontsize=18, numpoints=4)
    #["CTTHold", "CTTHnnaNT", "CTTHnnvNT", "CTTHnnm2NT", "CTTHnnmINT", "CTTHnnmINTnco2"]
    ax = fig.add_subplot(331)
    plot_one_subplot(aggregated_data, satz_step, pressure_plot,"CTTHold", "PPS-v2014")
    ax = fig.add_subplot(332)
    plot_one_subplot(aggregated_data, satz_step, pressure_plot, "CTTHnnaNT", "NN-AVHRR")
    ax = fig.add_subplot(333)
    plot_one_subplot(aggregated_data, satz_step, pressure_plot,"CTTHnna1NT", "NN-AVHRR1", legend=True)
    ax = fig.add_subplot(334)
    plot_one_subplot(aggregated_data, satz_step, pressure_plot,"CTTHnnvNT", "NN-VIIRS")
    ax = fig.add_subplot(335)
    plot_one_subplot(aggregated_data, satz_step, pressure_plot, "CTTHnnm2NT", "NN-MERSI-2")
    ax = fig.add_subplot(336)
    plot_one_subplot(aggregated_data, satz_step, pressure_plot,"CTTHnnmINT", "NN-MetImage")
    ax = fig.add_subplot(337)
    plot_one_subplot(aggregated_data, satz_step, pressure_plot,"CTTHnnmINTnco2", "NN-MetImage-NoCo2")
    ax = fig.add_subplot(338)
    plot_one_subplot(aggregated_data, satz_step, pressure_plot,"CTTH", "NN-AVHRR THR")
    ax = fig.add_subplot(339)
    plot_one_subplot(aggregated_data, satz_step, pressure_plot,"CTTHnna1", "NN-AVHRR1 THR")
    #ax.set_xlabel('pressure (hPa)')
    fig.text(0.5, 0.04, 'pressure (hPa)', ha='center',fontsize=20)
    fig.text(0.04, 0.5, 'percent of cloud top pressure results', va='center', rotation='vertical',fontsize=20)
    plt.savefig(PLOT_DIR + "ctth_satz_pressure_%s_m%s.png"%('modis',month), bbox_inches='tight')


    fig = plt.figure(figsize = (11,11))  
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')
    #plt.suptitle('Cloud top pressure dependence on satzenith angle', fontsize=12)
    def plot_one_subplot(aggregated_data, satz_step,pressure_plot,  ctth_label, desc, legend=False):
        #ax.text(420,11.5, desc, bbox={'facecolor':'w', 'edgecolor':'w','alpha':1.0, 'pad':6})
        plt.text(0.03, 0.90, desc, fontsize=18,transform=ax.transAxes, bbox=dict(facecolor='w', edgecolor='w', alpha=1.0))

        #ax.grid(True)
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        plot_00_20 = (aggregated_data[ctth_label][0] + 
                      aggregated_data[ctth_label][10])
        plot_20_40 = (aggregated_data[ctth_label][20] + 
                      aggregated_data[ctth_label][30])
        plot_40_60 = (aggregated_data[ctth_label][40] + 
                      aggregated_data[ctth_label][50])
        plot_60_80 = (aggregated_data[ctth_label][60] + 
                      aggregated_data[ctth_label][70])
        plot_calipso = aggregated_data["CALIPSO"][0]
        plot_00_20 = plot_00_20*100.0/np.sum(plot_00_20)
        plot_20_40 = plot_20_40*100.0/np.sum(plot_20_40)
        plot_40_60 = plot_40_60*100.0/np.sum(plot_40_60)
        plot_60_80 = plot_60_80*100.0/np.sum(plot_60_80)
        plot_calipso = plot_calipso*100.0/np.sum(plot_calipso) 
        plt.plot(pressure_plot[:-1], plot_calipso, 'b', alpha=0.2, label="CALIPSO") 
        plt.plot(pressure_plot[:-1], plot_00_20, '-k', label="satellite zenith angle: 0-20")
        plt.plot(pressure_plot[:-1], plot_20_40, '--k', label="satellite zenith angle: 20-40")
        plt.plot(pressure_plot[:-1], plot_40_60, '-.k', label="satellite zenith angle: 40-60")
        plt.plot(pressure_plot[:-1], plot_60_80, ':k', label="satellite zenith angle: 60-80")
        ax.fill(pressure_plot[:-1], plot_calipso, "b", alpha=0.2) 

        #ax.set_ylabel('percent of ctth results')
        ax.set_ylim(0,13)
        ax.set_xlim(70,1150)
        if legend:
            #my_legend = plt.legend(loc='upper right',bbox_to_anchor=(1.2, 1.2), framealpha=1.0,edgecolor='w')
            my_legend = plt.legend(loc='center left', bbox_to_anchor=(0.9, 0.6), framealpha=1.0,edgecolor='w',fontsize=18, numpoints=4)
            #frame = my_legend.get_frame()
            #frame.set_facecolor('white')
    #["CTTHold", "CTTHnnaNT", "CTTHnnvNT", "CTTHnnm2NT", "CTTHnnmINT", "CTTHnnmINTnco2"]

    ax = fig.add_subplot(331)
    plot_one_subplot(aggregated_data, satz_step, pressure_plot,"CTTHold", "(a) PPS-v2014")
    ax = fig.add_subplot(332)
    plot_one_subplot(aggregated_data, satz_step, pressure_plot,"MODIS-C6", "(b) MODIS-C6",legend=3)
    ax = fig.add_subplot(334)
    plot_one_subplot(aggregated_data, satz_step, pressure_plot, "CTTHnnaNT", "(c) NN-AVHRR")
    ax = fig.add_subplot(335)
    plot_one_subplot(aggregated_data, satz_step, pressure_plot,"CTTHnna1NT", "(d) NN-AVHRR1")
    ax = fig.add_subplot(336)
    plot_one_subplot(aggregated_data, satz_step, pressure_plot,"CTTHnnvNT", "(e) NN-VIIRS")
    ax = fig.add_subplot(337)
    plot_one_subplot(aggregated_data, satz_step, pressure_plot, "CTTHnnm2NT", "(f) NN-MERSI-2")
    ax = fig.add_subplot(338)
    plot_one_subplot(aggregated_data, satz_step, pressure_plot,"CTTHnnmINTnco2", "(g) NN-MetImage-$NoCO_2$")
    ax = fig.add_subplot(339)
    plot_one_subplot(aggregated_data, satz_step, pressure_plot,"CTTHnnmINT", "(h) NN-MetImage")

    #ax.set_xlabel('pressure (hPa)')
    fig.text(0.5, 0.04, 'Pressure (hPa)', ha='center',fontsize=18)
    fig.text(0.04, 0.5, 'Percent of cloud top pressure results', va='center', rotation='vertical',fontsize=18)
    plt.savefig(PLOT_DIR + "fig_01ctth_satz_pressure_%s_for_art_m%s.png"%('modis',month), bbox_inches='tight')
    plt.savefig(PLOT_DIR + "fig_01ctth_satz_pressure_%s_for_art_m%s.pdf"%('modis',month), bbox_inches='tight')
def investigate_nn_ctth_satz():
    month = '12'
    month_path = month
    month = '020406081012'
    month_path = '*'
    PLOT_DIR = "/home/a001865/PICTURES_FROM_PYTHON/CTTH_PLOTS/"
    out_filenames = glob("/home/a001865/DATA_MISC/satz_modis_01_investigation/*s_%s_*.txt"%(month_path))
    plot_pressure_modis_01(out_filenames,PLOT_DIR,month)
    out_filenames = glob("/home/a001865/git/atrain_match/satz_statistics_and_modis.txt")
    plot_pressure_local(out_filenames,PLOT_DIR)
if __name__ == "__main__":
    investigate_nn_ctth_satz()

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

# Author(s)
# Nina Hakansson

import numpy as np
from glob import glob
import os
try:
    from matchobject_io import (read_truth_imager_match_obj,
                                write_truth_imager_match_obj)
except:
    from atrain_match.matchobject_io import (read_truth_imager_match_obj, 
                                             write_truth_imager_match_obj)

import time
import pdb

instrument = "avhrr"
truth = "calipso"
version = "v2018"

SATELLITES = ["noaa18", "noaa19", "noaa17", "metopa", "metopb"]
YEAR_LIST = ['2006', '2007', '2008', '2009', '2010', '2011', '2012', '2013', '2014', '2015']
BASE_DIR = '/nobackup/smhid14/sm_erjoh/pps/Atrain_match/patchCMSAF_Oct2020_qualityflagcheck/Reshaped_Files'
# BASE_DIR = "/home/a001865/DATA_MISC/reshaped_files_from_kg_temp"
MAKE_EXTRA_CHECK = False

SETTINGS = {"WRITE_ONLY_THE_MOST_IMPORTANT_STUFF_TO_FILE": False}


def remove_doubles(mObj, mObj2):
    non_doubles = [ind for ind in range(len(mObj.calipso.sec_1970))
                   if mObj.calipso.sec_1970[ind] not in mObj2.calipso.sec_1970]

    if MAKE_EXTRA_CHECK:
        # Takes a lot of time ....
        doubles_time_and_id = [(mObj.calipso.sec_1970[ind], mObj.calipso.profile_id[ind, 0])
                               for ind in range(len(mObj.calipso.sec_1970)) if ind not in non_doubles]

        extra_check = [time_and_id for time_and_id in doubles_time_and_id if time_and_id in zip(
            mObj2.calipso.sec_1970, mObj2.calipso.profile_id[:, 0])]
        if len(extra_check) != len(mObj.calipso.sec_1970) - len(non_doubles):
            print("Some points identified as doubles does not have "
                  "the same time and profile id as ther corresponding double"
                  "There is something wrong (in the code) here!")
            raise ValueError

    print("Found {:d} doubles (out of {:d}) in the new object".format(
        len(mObj.calipso.sec_1970) - len(non_doubles), len(mObj.calipso.sec_1970)))
    
    #: All are doublets
    if len(non_doubles) == 0:
        return -1
    else:
        mObj = mObj.extract_elements(idx=np.array(non_doubles))
        return mObj
    # import pdb;pdb.set_trace()
        

if __name__ == '__main__':
    '''
    Creates monthly merged reshape files.
    BASE_DIR is the dir where original reshape files are stored
    Merged results will be in BASE_DIR + _merged_{truth}/{satellite}/
    '''
    tic = time.time()
    files_cant_read = []
    files_all_doubles = []
    match_calipso_merged = None
#     corrupt_files = ['5km_noaa18_20061101_2024_07479_calipso_avhrr_match.h5', '5km_noaa18_20080729_1251_16444_calipso_avhrr_match.h5', '5km_noaa18_20100709_1525_26460_calipso_avhrr_match.h5', '5km_noaa18_20150801_2020_52552_calipso_avhrr_match.h5', '5km_noaa19_20101102_0513_08937_calipso_avhrr_match.h5', '5km_noaa19_20111107_1930_14163_calipso_avhrr_match.h5', '5km_noaa19_20130103_0509_20122_calipso_avhrr_match.h5', '5km_noaa19_20130103_0656_20123_calipso_avhrr_match.h5', '5km_noaa19_20140129_0656_25639_calipso_avhrr_match.h5', '5km_noaa19_20150226_1002_31186_calipso_avhrr_match.h5', '5km_noaa19_20151117_0919_34911_calipso_avhrr_match.h5']
    corrupt_files = ['5km_noaa18_20061101_2024_07479_calipso_avhrr_match.h5', '5km_noaa18_20061112_1643_07631_calipso_avhrr_match.h5', '5km_noaa18_20061121_0705_07753_calipso_avhrr_match.h5', '5km_noaa18_20061125_1932_07816_calipso_avhrr_match.h5', '5km_noaa18_20061202_0513_07907_calipso_avhrr_match.h5', '5km_noaa18_20061202_0652_07908_calipso_avhrr_match.h5', '5km_noaa18_20061230_1023_08305_calipso_avhrr_match.h5', '5km_noaa18_20070127_1845_08704_calipso_avhrr_match.h5', '5km_noaa18_20070129_2013_08734_calipso_avhrr_match.h5', '5km_noaa18_20070216_1032_08982_calipso_avhrr_match.h5', '5km_noaa18_20070218_1333_09012_calipso_avhrr_match.h5', '5km_noaa18_20070314_1105_09349_calipso_avhrr_match.h5', '5km_noaa18_20070316_1226_09378_calipso_avhrr_match.h5', '5km_noaa18_20070318_2021_09411_calipso_avhrr_match.h5', '5km_noaa18_20070331_1941_09593_calipso_avhrr_match.h5', '5km_noaa18_20070407_0701_09685_calipso_avhrr_match.h5', '5km_noaa18_20070426_2021_09961_calipso_avhrr_match.h5', '5km_noaa18_20070501_0614_10023_calipso_avhrr_match.h5', '5km_noaa18_20070512_0223_10176_calipso_avhrr_match.h5', '5km_noaa18_20070525_0337_10360_calipso_avhrr_match.h5', '5km_noaa18_20070527_0829_10391_calipso_avhrr_match.h5', '5km_noaa18_20070527_1005_10392_calipso_avhrr_match.h5', '5km_noaa18_20070609_0944_10575_calipso_avhrr_match.h5', '5km_noaa18_20070701_0353_10882_calipso_avhrr_match.h5', '5km_noaa18_20070913_0248_11925_calipso_avhrr_match.h5', '5km_noaa18_20071121_1920_12907_calipso_avhrr_match.h5', '5km_noaa18_20071204_2034_13092_calipso_avhrr_match.h5', '5km_noaa18_20071209_0621_13154_calipso_avhrr_match.h5', '5km_noaa18_20080104_0841_13522_calipso_avhrr_match.h5', '5km_noaa18_20080106_1317_13553_calipso_avhrr_match.h5', '5km_noaa18_20080108_1605_13583_calipso_avhrr_match.h5', '5km_noaa18_20080221_0515_14197_calipso_avhrr_match.h5', '5km_noaa18_20080223_0959_14228_calipso_avhrr_match.h5', '5km_noaa18_20080409_0517_14874_calipso_avhrr_match.h5', '5km_noaa18_20080422_0623_15058_calipso_avhrr_match.h5', '5km_noaa18_20080428_2022_15151_calipso_avhrr_match.h5', '5km_noaa18_20080503_0609_15213_calipso_avhrr_match.h5', '5km_noaa18_20080518_1159_15428_calipso_avhrr_match.h5', '5km_noaa18_20080613_1715_15797_calipso_avhrr_match.h5', '5km_noaa18_20080729_1251_16444_calipso_avhrr_match.h5', '5km_noaa18_20080809_0738_16596_calipso_avhrr_match.h5', '5km_noaa18_20081113_0739_17950_calipso_avhrr_match.h5', '5km_noaa18_20081113_0915_17951_calipso_avhrr_match.h5', '5km_noaa18_20081207_1006_18290_calipso_avhrr_match.h5', '5km_noaa18_20081224_2025_18536_calipso_avhrr_match.h5', '5km_noaa18_20090124_1005_18967_calipso_avhrr_match.h5', '5km_noaa18_20090126_1306_18997_calipso_avhrr_match.h5', '5km_noaa18_20090410_2007_20045_calipso_avhrr_match.h5', '5km_noaa18_20090426_0538_20262_calipso_avhrr_match.h5', '5km_noaa18_20090509_0649_20446_calipso_avhrr_match.h5', '5km_noaa18_20090515_2041_20539_calipso_avhrr_match.h5', '5km_noaa18_20090522_0754_20630_calipso_avhrr_match.h5', '5km_noaa18_20090522_1111_20632_calipso_avhrr_match.h5', '5km_noaa18_20090615_0843_20969_calipso_avhrr_match.h5', '5km_noaa18_20090615_1019_20970_calipso_avhrr_match.h5', '5km_noaa18_20090615_1200_20971_calipso_avhrr_match.h5', '5km_noaa18_20090729_0221_21586_calipso_avhrr_match.h5', '5km_noaa18_20090906_0920_22140_calipso_avhrr_match.h5', '5km_noaa18_20091028_1959_22880_calipso_avhrr_match.h5', '5km_noaa18_20091108_1943_23035_calipso_avhrr_match.h5', '5km_noaa18_20091115_0837_23127_calipso_avhrr_match.h5', '5km_noaa18_20091121_2050_23219_calipso_avhrr_match.h5', '5km_noaa18_20091224_1959_23684_calipso_avhrr_match.h5', '5km_noaa18_20100109_0156_23899_calipso_avhrr_match.h5', '5km_noaa18_20100111_0319_23928_calipso_avhrr_match.h5', '5km_noaa18_20100111_0655_23930_calipso_avhrr_match.h5', '5km_noaa18_20100228_0831_24608_calipso_avhrr_match.h5', '5km_noaa18_20100317_2024_24855_calipso_avhrr_match.h5', '5km_noaa18_20100424_0144_25380_calipso_avhrr_match.h5', '5km_noaa18_20100426_0648_25411_calipso_avhrr_match.h5', '5km_noaa18_20100709_1525_26460_calipso_avhrr_match.h5', '5km_noaa18_20100727_0528_26708_calipso_avhrr_match.h5', '5km_noaa18_20100822_0936_27077_calipso_avhrr_match.h5', '5km_noaa18_20100902_0912_27232_calipso_avhrr_match.h5', '5km_noaa18_20101109_0701_28190_calipso_avhrr_match.h5', '5km_noaa18_20101113_1925_28253_calipso_avhrr_match.h5', '5km_noaa18_20101207_1828_28591_calipso_avhrr_match.h5', '5km_noaa18_20101207_2017_28593_calipso_avhrr_match.h5', '5km_noaa18_20101229_1936_28902_calipso_avhrr_match.h5', '5km_noaa18_20110122_2028_29242_calipso_avhrr_match.h5', '5km_noaa18_20110202_2005_29396_calipso_avhrr_match.h5', '5km_noaa18_20110209_1221_29491_calipso_avhrr_match.h5', '5km_noaa18_20110503_2051_30667_calipso_avhrr_match.h5', '5km_noaa18_20110514_2029_30821_calipso_avhrr_match.h5', '5km_noaa18_20111207_0510_33733_calipso_avhrr_match.h5', '5km_noaa18_20111220_1257_33921_calipso_avhrr_match.h5', '5km_noaa18_20120102_1406_34105_calipso_avhrr_match.h5', '5km_noaa18_20120611_0817_36373_calipso_avhrr_match.h5', '5km_noaa18_20120816_1133_37306_calipso_avhrr_match.h5', '5km_noaa18_20121028_0757_38334_calipso_avhrr_match.h5', '5km_noaa18_20121110_1402_38521_calipso_avhrr_match.h5', '5km_noaa18_20130102_1413_39269_calipso_avhrr_match.h5', '5km_noaa18_20130305_1111_40142_calipso_avhrr_match.h5', '5km_noaa18_20130423_0334_40829_calipso_avhrr_match.h5', '5km_noaa18_20130716_0450_42015_calipso_avhrr_match.h5', '5km_noaa18_20130720_1112_42075_calipso_avhrr_match.h5', '5km_noaa18_20130809_0909_42356_calipso_avhrr_match.h5', '5km_noaa18_20130818_0720_42482_calipso_avhrr_match.h5', '5km_noaa18_20131222_1221_44263_calipso_avhrr_match.h5', '5km_noaa18_20140102_1339_44419_calipso_avhrr_match.h5', '5km_noaa18_20140131_0958_44826_calipso_avhrr_match.h5', '5km_noaa18_20140131_1316_44828_calipso_avhrr_match.h5', '5km_noaa18_20141219_1614_49374_calipso_avhrr_match.h5', '5km_noaa18_20141226_0814_49468_calipso_avhrr_match.h5', '5km_noaa18_20150202_0736_50004_calipso_avhrr_match.h5', '5km_noaa18_20150502_0911_51261_calipso_avhrr_match.h5', '5km_noaa18_20150801_2020_52552_calipso_avhrr_match.h5', '5km_noaa19_20090316_2012_00541_calipso_avhrr_match.h5', '5km_noaa19_20090620_0722_01887_calipso_avhrr_match.h5', '5km_noaa19_20090701_0704_02042_calipso_avhrr_match.h5', '5km_noaa19_20090714_0442_02224_calipso_avhrr_match.h5', '5km_noaa19_20090824_1118_02806_calipso_avhrr_match.h5', '5km_noaa19_20090828_1705_02865_calipso_avhrr_match.h5', '5km_noaa19_20091030_1945_03755_calipso_avhrr_match.h5', '5km_noaa19_20091119_0626_04030_calipso_avhrr_match.h5', '5km_noaa19_20091204_1359_04246_calipso_avhrr_match.h5', '5km_noaa19_20091204_1539_04247_calipso_avhrr_match.h5', '5km_noaa19_20091208_1946_04305_calipso_avhrr_match.h5', '5km_noaa19_20091215_0707_04397_calipso_avhrr_match.h5', '5km_noaa19_20091215_0849_04398_calipso_avhrr_match.h5', '5km_noaa19_20091215_1024_04399_calipso_avhrr_match.h5', '5km_noaa19_20091219_1933_04460_calipso_avhrr_match.h5', '5km_noaa19_20100104_0137_04676_calipso_avhrr_match.h5', '5km_noaa19_20100121_1043_04921_calipso_avhrr_match.h5', '5km_noaa19_20100123_1347_04951_calipso_avhrr_match.h5', '5km_noaa19_20100127_1934_05010_calipso_avhrr_match.h5', '5km_noaa19_20100203_1013_05104_calipso_avhrr_match.h5', '5km_noaa19_20100220_2036_05350_calipso_avhrr_match.h5', '5km_noaa19_20100303_1348_05501_calipso_avhrr_match.h5', '5km_noaa19_20100314_1013_05654_calipso_avhrr_match.h5', '5km_noaa19_20100725_0354_07526_calipso_avhrr_match.h5', '5km_noaa19_20100822_0919_07924_calipso_avhrr_match.h5', '5km_noaa19_20101102_0513_08937_calipso_avhrr_match.h5', '5km_noaa19_20101108_1905_09029_calipso_avhrr_match.h5', '5km_noaa19_20101110_2033_09059_calipso_avhrr_match.h5', '5km_noaa19_20101204_1945_09397_calipso_avhrr_match.h5', '5km_noaa19_20101211_0840_09489_calipso_avhrr_match.h5', '5km_noaa19_20110106_1545_09860_calipso_avhrr_match.h5', '5km_noaa19_20110115_0550_09981_calipso_avhrr_match.h5', '5km_noaa19_20110130_1321_10197_calipso_avhrr_match.h5', '5km_noaa19_20110212_1246_10380_calipso_avhrr_match.h5', '5km_noaa19_20110319_1002_10872_calipso_avhrr_match.h5', '5km_noaa19_20110405_2023_11118_calipso_avhrr_match.h5', '5km_noaa19_20110504_0521_11518_calipso_avhrr_match.h5', '5km_noaa19_20110523_2026_11795_calipso_avhrr_match.h5', '5km_noaa19_20110616_1930_12132_calipso_avhrr_match.h5', '5km_noaa19_20111014_2027_13826_calipso_avhrr_match.h5', '5km_noaa19_20111021_0740_13917_calipso_avhrr_match.h5', '5km_noaa19_20111107_1930_14163_calipso_avhrr_match.h5', '5km_noaa19_20111230_0707_14904_calipso_avhrr_match.h5', '5km_noaa19_20120101_1004_14934_calipso_avhrr_match.h5', '5km_noaa19_20120103_1612_14966_calipso_avhrr_match.h5', '5km_noaa19_20120209_1951_15490_calipso_avhrr_match.h5', '5km_noaa19_20120218_1000_15611_calipso_avhrr_match.h5', '5km_noaa19_20120402_0532_16229_calipso_avhrr_match.h5', '5km_noaa19_20120415_0642_16413_calipso_avhrr_match.h5', '5km_noaa19_20120430_1709_16630_calipso_avhrr_match.h5', '5km_noaa19_20120613_0801_17246_calipso_avhrr_match.h5', '5km_noaa19_20120613_0937_17247_calipso_avhrr_match.h5', '5km_noaa19_20120613_1118_17248_calipso_avhrr_match.h5', '5km_noaa19_20120707_1024_17586_calipso_avhrr_match.h5', '5km_noaa19_20121102_1924_19255_calipso_avhrr_match.h5', '5km_noaa19_20121107_0516_19318_calipso_avhrr_match.h5', '5km_noaa19_20121124_1402_19563_calipso_avhrr_match.h5', '5km_noaa19_20121124_1542_19564_calipso_avhrr_match.h5', '5km_noaa19_20121207_1958_19750_calipso_avhrr_match.h5', '5km_noaa19_20121218_1941_19905_calipso_avhrr_match.h5', '5km_noaa19_20130103_0509_20122_calipso_avhrr_match.h5', '5km_noaa19_20130103_0656_20123_calipso_avhrr_match.h5', '5km_noaa19_20130111_2208_20245_calipso_avhrr_match.h5', '5km_noaa19_20130122_2008_20399_calipso_avhrr_match.h5', '5km_noaa19_20130127_0741_20462_calipso_avhrr_match.h5', '5km_noaa19_20130127_0918_20463_calipso_avhrr_match.h5', '5km_noaa19_20130301_0648_20927_calipso_avhrr_match.h5', '5km_noaa19_20130307_2039_21020_calipso_avhrr_match.h5', '5km_noaa19_20130318_2021_21175_calipso_avhrr_match.h5', '5km_noaa19_20130728_1150_23032_calipso_avhrr_match.h5', '5km_noaa19_20130819_1113_23342_calipso_avhrr_match.h5', '5km_noaa19_20131011_0826_24088_calipso_avhrr_match.h5', '5km_noaa19_20131102_0748_24398_calipso_avhrr_match.h5', '5km_noaa19_20131119_1932_24644_calipso_avhrr_match.h5', '5km_noaa19_20131124_0523_24707_calipso_avhrr_match.h5', '5km_noaa19_20131124_0710_24708_calipso_avhrr_match.h5', '5km_noaa19_20131126_1328_24740_calipso_avhrr_match.h5', '5km_noaa19_20131203_0337_24833_calipso_avhrr_match.h5', '5km_noaa19_20131205_0504_24862_calipso_avhrr_match.h5', '5km_noaa19_20131211_2041_24956_calipso_avhrr_match.h5', '5km_noaa19_20131218_1109_25049_calipso_avhrr_match.h5', '5km_noaa19_20131229_1047_25204_calipso_avhrr_match.h5', '5km_noaa19_20131229_1231_25205_calipso_avhrr_match.h5', '5km_noaa19_20140105_0238_25298_calipso_avhrr_match.h5', '5km_noaa19_20140122_1632_25546_calipso_avhrr_match.h5', '5km_noaa19_20140127_0154_25608_calipso_avhrr_match.h5', '5km_noaa19_20140129_0656_25639_calipso_avhrr_match.h5', '5km_noaa19_20140131_1313_25671_calipso_avhrr_match.h5', '5km_noaa19_20140209_0630_25794_calipso_avhrr_match.h5', '5km_noaa19_20140325_1157_26418_calipso_avhrr_match.h5', '5km_noaa19_20140501_2004_26944_calipso_avhrr_match.h5', '5km_noaa19_20140528_0700_27318_calipso_avhrr_match.h5', '5km_noaa19_20140729_0534_28192_calipso_avhrr_match.h5', '5km_noaa19_20140731_0841_28222_calipso_avhrr_match.h5', '5km_noaa19_20141012_0503_29250_calipso_avhrr_match.h5', '5km_noaa19_20141226_0949_30311_calipso_avhrr_match.h5', '5km_noaa19_20141226_1311_30313_calipso_avhrr_match.h5', '5km_noaa19_20141228_1734_30343_calipso_avhrr_match.h5', '5km_noaa19_20150108_1548_30498_calipso_avhrr_match.h5', '5km_noaa19_20150119_2020_30656_calipso_avhrr_match.h5', '5km_noaa19_20150208_1951_30937_calipso_avhrr_match.h5', '5km_noaa19_20150226_1002_31186_calipso_avhrr_match.h5', '5km_noaa19_20150307_0638_31311_calipso_avhrr_match.h5', '5km_noaa19_20150311_1906_31374_calipso_avhrr_match.h5', '5km_noaa19_20150322_2032_31531_calipso_avhrr_match.h5', '5km_noaa19_20150331_1704_31655_calipso_avhrr_match.h5', '5km_noaa19_20150707_0349_33031_calipso_avhrr_match.h5', '5km_noaa19_20150707_0726_33033_calipso_avhrr_match.h5', '5km_noaa19_20151026_0448_34598_calipso_avhrr_match.h5', '5km_noaa19_20151028_0802_34628_calipso_avhrr_match.h5', '5km_noaa19_20151106_0614_34754_calipso_avhrr_match.h5', '5km_noaa19_20151117_0919_34911_calipso_avhrr_match.h5', '5km_noaa19_20151225_0515_35445_calipso_avhrr_match.h5']
    for satellite in SATELLITES:
        ROOT_DIR = BASE_DIR + "/{satellite}/5km/%s/%s/*%s%s*_*{truth}*.h5".format(
            satellite=satellite, truth=truth)
        print(ROOT_DIR)
        OUT_DIR_TEMPLATE = BASE_DIR + "_merged_{truth}/{satellite}/".format(
            truth=truth, satellite=satellite)
        outfile_template = "5km_{satellite}_%s%s01_0000_99999_{truth}_{instrument}_match.h5".format(
            satellite=satellite, truth=truth, instrument=instrument)
    
        for year in YEAR_LIST:
            for month in ["01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12"]:
                num_n = 0
                files_all_doubles_month = []
                files = sorted(glob(ROOT_DIR % (year, month,
                                                year, month)))
                if len(files) == 0:
                    continue
                print('\n Year = %s, mon = %s \n' %(year, month))
                OUT_DIR = OUT_DIR_TEMPLATE
                if not os.path.exists(OUT_DIR):
                    os.makedirs(OUT_DIR)
                for filename in files:
                    if os.path.basename(filename) in corrupt_files:
                        print("File is corrupt. Do not use")
                        continue
                    print(os.path.basename(filename))
                    # match_calipso_new=read_truth_imager_match_obj(filename, truth=truth)
                    try:
                        match_calipso_new = read_truth_imager_match_obj(filename, truth=truth)
                    except:
                        print("problem with", os.path.basename(filename))
                        files_cant_read.append(filename)
                        continue
                        # if (match_calipso_new.cloudsat.RVOD_CWC_status is None or
                        # len(match_calipso_new.cloudsat.RVOD_CWC_status) != len(match_calipso_new.avhrr.cpp_lwp)):
                        #   print("Missing RVOD_CWC_status")
                        #   continue
                        
                    print("reading", os.path.basename(filename))
                    if match_calipso_merged is None:
                        match_calipso_merged = match_calipso_new
                    else:
                        match_calipso_new = remove_doubles(match_calipso_new, match_calipso_merged)
                        if match_calipso_new == -1:
                            files_all_doubles.append(filename)
                            files_all_doubles_month.append(filename)
                            continue
                        match_calipso_merged = match_calipso_merged + match_calipso_new
                    num_n += 1
    
                if num_n > 0:
                    filename_merged = outfile_template % (year, month)
                    outfile = os.path.join(OUT_DIR, filename_merged)
    
                    print(len(match_calipso_merged.calipso.sec_1970))
                    write_truth_imager_match_obj(
                        outfile, match_calipso_merged, SETTINGS, imager_obj_name='pps')
                    match_calipso_merged = None
                    print("File with all doublets %s files" %len(files_all_doubles_month))
    
    print("Could not read %s files" %len(files_cant_read))
    print("File with all doublets %s files, total" %len(files_all_doubles))
    print(time.time()-tic)
    if (len(files_cant_read) != 0) or (len(files_all_doubles) != 0):
        pdb.set_trace()

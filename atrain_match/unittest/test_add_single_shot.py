# -*- coding: utf-8 -*-
#
# test_match.py
"""for the match module match.py """
#
import numpy as np
import unittest
from truths.calipso import  (addSingleShotTo5km)
from matchobject_io import CalipsoObject


def addSingleShotTo5km_old(Obj5): # Valid only for CALIPSO-CALIOP version 4.10
    # weakness or bug in the CALIPSO retrieval of clouds below 4 km
    for i in range(Obj5.profile_utc_time.shape[0]):
        cfc = Obj5.Number_cloudy_single_shots[i]/15.0 
        if cfc == 1.0:
            cfc = 0.99
        if (Obj5.number_layers_found[i] > 0):
            Obj5.cloud_fraction[i] = 1.0
        if ((cfc > 0.01) and (Obj5.number_layers_found[i] == 0)):
            Obj5.number_layers_found[i] = 1
            Obj5.layer_top_altitude[i, 0] = Obj5.Average_cloud_top_single_shots[i]
            Obj5.layer_base_altitude[i, 0] = Obj5.Average_cloud_base_single_shots[i]
            Obj5.feature_optical_depth_532[i, 0] = 1.0
            Obj5.cloud_fraction[i] = cfc
    return Obj5


class test_addSingleShot(unittest.TestCase): 

    def setUp(self):
        
        self.obt5 = CalipsoObject()
        self.obt5.profile_utc_time = np.zeros((9,3)) 
        self.obt5.layer_top_altitude = np.zeros((9,10)) -9
        self.obt5.layer_base_altitude = np.zeros((9,10)) -9
        self.obt5.layer_top_pressure = np.zeros((9,10)) -9
        self.obt5.layer_base_pressure = np.zeros((9,10)) -9
        self.obt5.feature_optical_depth_532 = np.zeros((9,10)) -9
        self.obt5.number_layers_found = np.array([0,0,0,0,0,0,0,3,3])
        self.obt5.cloud_fraction = np.array([0,0,0,0,0,0,0,1,1])
        self.obt5.feature_classification_flags = np.zeros((9,10)) -9
        self.obt5.Number_cloudy_single_shots =  np.array([15, 0, 1,  10, 10,   0, 5,  14, 0]).ravel()
        self.obt5.Average_cloud_top_single_shots = np.array([15000, 0, 1000,  1000, 8000,   0, 500,  1400, 0]).ravel()
        self.obt5.Average_cloud_base_single_shots = np.array([14000, 0, 500,  50, 3000,   0, 50,  1400, 0]).ravel()
        self.obt5.layer_top_altitude[:,0] =        np.array([-9, -9, -9, -9, -9,   -9, -9,  2.2, 5.0]).ravel()
        self.obt5.layer_base_altitude[:,0] =       np.array([-9, -9, -9, -9, -9,   -9, -9,  2.0, 3.1]).ravel()
        self.obt5.feature_optical_depth_532[:,0] = np.array([-9, -9, -9,  -9, -9,  -9, -9, 10.0, 0.3]).ravel()
        self.obt5.feature_classification_flags[:,0] = np.array([-9,   -9, -9, -9,   -9,   -9, -9, 2, 3]).ravel()


    def test_new_code(self):
        out1 = addSingleShotTo5km(self.obt5)
        out2 = addSingleShotTo5km_old(self.obt5)
        #This is what I think we should do

        self.assertTrue((np.equal(out1.layer_top_altitude,out2.layer_top_altitude)).all())
        self.assertTrue((np.equal(out1.number_layers_found,out2.number_layers_found)).all())
        self.assertTrue((np.equal(out1.layer_base_altitude,out2.layer_base_altitude)).all())
        self.assertTrue((np.equal(out1.feature_optical_depth_532,out2.feature_optical_depth_532)).all())
        self.assertTrue((np.equal(out1.cloud_fraction,out2.cloud_fraction)).all())


def suite():
    """The suite for test_utils.
    """
    loader = unittest.TestLoader()
    mysuite = unittest.TestSuite()
    mysuite.addTest(loader.loadTestsFromTestCase(test_addSingleShot))
    return mysuite

if __name__ == "__main__":
    unittest.TextTestRunner(verbosity=2).run(suite())

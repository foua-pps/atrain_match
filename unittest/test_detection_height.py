# -*- coding: utf-8 -*-
#
# test_match.py
"""for the match module match.py """
#
import numpy as np
import unittest
from calipso import  (detection_height_from_5km_data)
from matchobject_io import CalipsoObject
class test_detection_height(unittest.TestCase): 

    def setUp(self):
        
        self.obt5 = CalipsoObject()
        self.obt1 = CalipsoObject()
        self.obt5.profile_utc_time = np.zeros((7,3)) 
        self.obt1.profile_utc_time = np.zeros((35,1))-10
        self.obt1.profile_utc_time[0:7] = 0
        self.obt1.number_layers_found = np.ones((35,1))
        self.obt5.layer_top_altitude = np.zeros((7,10)) -9999
        self.obt5.layer_base_altitude = np.zeros((7,10)) -9999
        self.obt5.feature_optical_depth_532 = np.zeros((7,10)) -9999
        self.obt5.number_layers_found = np.array([1,1,1,3,3,3,3])
        self.obt1.total_optical_depth_5km = np.repeat(np.array([0.1, 9.0 , 0.3, 0.2, 3.6, 10.4, 0.9]), 5, axis=0)

        self.obt5.layer_top_altitude[:,0] = np.array([-9999,   9.2,    8.3,   7.6, 5.3, 2.2, 5.0]).ravel()
        self.obt5.layer_top_altitude[:,1] = np.array([-9999,  -9999, -9999,   6.1, 4.3, 1.2, 4.9]).ravel()
        self.obt5.layer_top_altitude[:,2] = np.array([-9999,  -9999, -9999,   5.6, 3.3, 0.6, 4.8]).ravel()
        self.obt5.layer_base_altitude[:,0] = np.array([-9999,   8.2,   7.3,   6.6, 1.3, 2.0, 3.1]).ravel()
        self.obt5.layer_base_altitude[:,1] = np.array([-9999, -9999, -9999,   4.6, 1.2, 0.2, 2.1]).ravel()
        self.obt5.layer_base_altitude[:,2] = np.array([-9999, -9999, -9999, -9999, 1.1, 0.1, 1.1]).ravel()

        self.obt5.feature_optical_depth_532[:,0] = np.array([-9999,     8.2, 0.3,   0.1, 1.3, 10.0, 0.3]).ravel()
        self.obt5.feature_optical_depth_532[:,1] = np.array([-9999, -9999, -9999,   0.1, 1.2,  0.2, 0.3]).ravel()
        self.obt5.feature_optical_depth_532[:,2] = np.array([-9999, -9999, -9999, -9999, 1.1,  0.1, 0.3]).ravel()



    def test_detection(self):
        calipso = detection_height_from_5km_data(self.obt1, self.obt5,
                                                     1,33)
        #print calipso.detection_height_5km                                                    
        self.assertEqual(calipso.detection_height_5km[0], -9)
        self.assertTrue(np.abs(calipso.detection_height_5km[4]-9.078)<0.01)
        self.assertTrue(np.abs(calipso.detection_height_5km[10]-7.3)<0.01)
        self.assertTrue(np.abs(calipso.detection_height_5km[15]-0.0)<0.01)
        self.assertTrue(np.abs(calipso.detection_height_5km[20]-3.353)<0.01)
        self.assertTrue(np.abs(calipso.detection_height_5km[25]-2.18)<0.01)
        self.assertTrue(np.abs(calipso.detection_height_5km[30]-0.0)<0.01)

def suite():
    """The suite for test_utils.
    """
    loader = unittest.TestLoader()
    mysuite = unittest.TestSuite()
    mysuite.addTest(loader.loadTestsFromTestCase(test_detection_height))
    return mysuite

if __name__ == "__main__":
    unittest.TextTestRunner(verbosity=2).run(suite())

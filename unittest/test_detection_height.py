# -*- coding: utf-8 -*-
#
# test_match.py
"""for the match module match.py """
#
import numpy as np
import unittest
from cloudsat_calipso_avhrr_prepare import  (detection_height_from_5km_data,
                                             CalipsoCloudOpticalDepth_new)
from matchobject_io import CalipsoObject


def CalipsoCloudOpticalDepth(cloud_top, cloud_base, optical_depth, cloud_fraction, fcf, min_optical_depth):
    from config import MIN_OPTICAL_DEPTH
    new_cloud_top = np.ones(cloud_top.shape,'d')*np.min(cloud_top)
    new_cloud_base = np.ones(cloud_base.shape,'d')*np.min(cloud_base)
    new_cloud_fraction = np.zeros(cloud_fraction.shape,'d')
    new_fcf = np.ones(fcf.shape).astype(fcf.dtype)
    for pixel_i in range(optical_depth.shape[0]):
        depthsum = 0 #Used to sum the optical_depth
        for layer_j in range(optical_depth.shape[1]):
            # Just stops the for loop when there are no more valid value 
            if optical_depth[pixel_i, layer_j] < 0:
                break
            else:
                depthsum = depthsum + optical_depth[pixel_i, layer_j]
                if depthsum >= min_optical_depth:
                    new_cloud_top[pixel_i, 0:(optical_depth.shape[1]-layer_j)] = cloud_top[pixel_i, layer_j:]
                    new_cloud_base[pixel_i, 0:(optical_depth.shape[1]-layer_j)] = cloud_base[pixel_i, layer_j:]
                    new_cloud_fraction[pixel_i] = cloud_fraction[pixel_i] 
                    new_fcf[pixel_i, 0:(fcf.shape[1]-layer_j)] = fcf[pixel_i, layer_j:]
                    for k in range(new_cloud_top.shape[1]):
                        if new_cloud_top[pixel_i, k] < 0:
                            break
                        else:
                            new_cloud_top[pixel_i, k] = new_cloud_base[pixel_i, k] + \
                            ((new_cloud_top[pixel_i, k] - new_cloud_base[pixel_i, k]) * 1/2.)
                    break
    new_validation_height = new_cloud_top[:,0].copy()
    new_validation_height[new_validation_height>=0] *= 1000
    new_validation_height[new_validation_height<0] = -9
    return new_cloud_top, new_cloud_base, new_cloud_fraction, new_fcf, new_validation_height




class test_detection_height(unittest.TestCase): 

    def setUp(self):
        
        self.obt5 = CalipsoObject()
        self.obt1 = CalipsoObject()
        self.obt5.profile_utc_time = np.zeros((7,3)) 
        self.obt1.profile_utc_time = np.zeros((35,1))-10
        self.obt1.profile_utc_time[0:7] = 0
        self.obt1.number_layers_found = np.ones((35,1))
        self.obt5.layer_top_altitude = np.zeros((7,10)) -9
        self.obt5.layer_base_altitude = np.zeros((7,10)) -9
        self.obt5.layer_top_pressure = np.zeros((7,10)) -9
        self.obt5.layer_base_pressure = np.zeros((7,10)) -9
        self.obt5.feature_optical_depth_532 = np.zeros((7,10)) -9
        self.obt5.number_layers_found = np.array([1,1,1,3,3,3,3])
        self.obt5.cloud_fraction = np.array([1,1,1,1,1,1,1])
        self.obt5.feature_classification_flags = np.zeros((7,10)) -9
        self.obt1.total_optical_depth_5km = np.repeat(np.array([0.1, 9.0 , 0.3, 0.2, 3.6, 10.4, 0.9]), 5, axis=0)

        self.obt5.layer_top_altitude[:,0] =        np.array([-9, 9.2, 8.3,   7.6, 5.3,  2.2, 5.0]).ravel()
        self.obt5.layer_top_altitude[:,1] =        np.array([-9,  -9,  -9,   6.1, 4.3,  1.2, 4.9]).ravel()
        self.obt5.layer_top_altitude[:,2] =        np.array([-9,  -9,  -9,   5.6, 3.3,  0.6, 4.8]).ravel()
        self.obt5.layer_base_altitude[:,0] =       np.array([-9, 8.2, 7.3,   6.6, 1.3,  2.0, 3.1]).ravel()
        self.obt5.layer_base_altitude[:,1] =       np.array([-9,  -9,  -9,   4.6, 1.2,  0.2, 2.1]).ravel()
        self.obt5.layer_base_altitude[:,2] =       np.array([-9,  -9,  -9,    -9, 1.1 , 0.1, 1.1]).ravel()

        self.obt5.feature_optical_depth_532[:,0] = np.array([-9, 8.2, 0.3,   0.1, 1.3, 10.0, 0.3]).ravel()
        self.obt5.feature_optical_depth_532[:,1] = np.array([-9,  -9,  -9,   0.1, 1.2,  0.2, 0.3]).ravel()
        self.obt5.feature_optical_depth_532[:,2] = np.array([-9,  -9,  -9,    -9, 1.1,  0.1, 0.3]).ravel()

        self.obt5.feature_classification_flags[:,0] = np.array([-9,   8,   7,   6, 1, 2, 3]).ravel()
        self.obt5.feature_classification_flags[:,1] = np.array([-9, -9, -9,   4, 1, 2, 2]).ravel()
        self.obt5.feature_classification_flags[:,2] = np.array([-9, -9, -9, -9, 9, 1, 1]).ravel()




    def test_detection(self):
        calipso = detection_height_from_5km_data(self.obt1, self.obt5, limit_ctop=1.0)
        print "dh", calipso.detection_height_5km[0::5]                                                    
        self.assertEqual(calipso.detection_height_5km[0], -9)
        self.assertTrue(np.abs(calipso.detection_height_5km[5]-9.08*1000)<100)
        self.assertTrue(np.abs(calipso.detection_height_5km[10]-7.3*1000)<100)
        self.assertTrue(np.abs(calipso.detection_height_5km[15]-6.6*1000)<100)
        self.assertTrue(np.abs(calipso.detection_height_5km[20]-2.22*1000)<100)
        self.assertTrue(np.abs(calipso.detection_height_5km[25]-2.18*1000)<100)
        self.assertTrue(np.abs(calipso.detection_height_5km[30]-3.1*1000)<100)
    def test_ninas_code(self):
        out1 = CalipsoCloudOpticalDepth_new(self.obt5, 0.5, use_old_method=True,
                                             limit_ctop=1.0)
        out2 = CalipsoCloudOpticalDepth(self.obt5.layer_top_altitude, 
                                        self.obt5.layer_base_altitude, 
                                        self.obt5.feature_optical_depth_532,
                                        self.obt5.cloud_fraction, 
                                        self.obt5.feature_classification_flags, 
                                        0.5) 
        #This is what I think we should do
        out3 = CalipsoCloudOpticalDepth_new(self.obt5, 0.5, use_old_method=False,
                                             limit_ctop=0.1)
        self.assertTrue((np.equal(out1[0],out2[0])).all())
        self.assertTrue((np.equal(out1[2],out2[2])).all())
        self.assertTrue((np.equal(out1[3],out2[3])).all())
        self.assertTrue((np.equal(out1[4],out2[4])).all())
        self.assertTrue((np.equal(out1[1],out2[1])).all())


def suite():
    """The suite for test_utils.
    """
    loader = unittest.TestLoader()
    mysuite = unittest.TestSuite()
    mysuite.addTest(loader.loadTestsFromTestCase(test_detection_height))
    return mysuite

if __name__ == "__main__":
    unittest.TextTestRunner(verbosity=2).run(suite())

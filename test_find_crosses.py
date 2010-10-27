import unittest
import find_crosses
from datetime import datetime

TEST_CROSS = find_crosses.Cross('NOAA18', 'CALIPSO',
                                     datetime(2010, 1, 7, 0, 59, 11),
                                     datetime(2010, 1, 7, 0,  51, 58),
                                     25.21, 54.17)
HRPT_FILENAME = "%s/201001/%s" % (find_crosses.HRPT_DIR, 
                                  'hrpt16_NOAA-18_07-JAN-2010_00:51:20.078_23870')
AVHRR_FILENAME = "%s/%s" % (find_crosses.SM_AVHRR_DIR,
                            'noaa18_20100107_0051_23871_satproj_00000_04989_avhrr.h5')


class Test_find_crosses(unittest.TestCase):
    
    
    def test_parse_line(self):
        correct_line = "  2010  1  7  0 59  11.4     2010  1  7  0 51  58.8        54.17    25.21"
        (time1, time2, lon, lat) = find_crosses.parse_line(correct_line)
        self.assertEquals(time1, datetime(2010, 1, 7, 0, 59, 11))
        self.assertEquals(time2, datetime(2010, 1, 7, 0,  51, 58))
        self.assertEquals(lon, 25.21)
        self.assertEquals(lat, 54.17)
    
    
    def test_parse_hrpt_filename(self):
        parsed = find_crosses.parse_hrpt_filename(HRPT_FILENAME)
        self.assertEquals(parsed['satellite'], 'NOAA-18')
        self.assertEquals(parsed['datetime'], datetime(2010, 1, 7, 0, 51, 20))
        self.assertEquals(parsed['orbit'], 23870)
    
    
    def test_find_hrpt_file(self):
        file_found = find_crosses.find_hrpt_file(TEST_CROSS)
        self.assertEquals(file_found, HRPT_FILENAME)
    
    
    def test_find_remapped_avhrr_file(self):
        file_found = find_crosses.find_remapped_avhrr_file(TEST_CROSS)
        self.assertEquals(file_found, AVHRR_FILENAME)


if __name__ == '__main__':
    unittest.main()


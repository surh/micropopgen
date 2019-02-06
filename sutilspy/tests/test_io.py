import unittest
import os
from sutilspy import io
#import shutils

class TestWriteTable(unittest.TestCase):
    """Test write_table"""
        
    def setUp(self):
        self.path = 'temp_test_sutilspy.io.write_table.txt'
        if os.path.isfile(self.path):
            raise FileExistsError("File {} already exists. Cannot overwrite".format(self.path))
    
    def tearDown(self):
        os.unlink(self.path)
    
    def test_write_table(self):
        """Test correct number of lines"""
        
        rows = ['row1','row2']
        header = 'header'
        delimiter = ','
        verbose = True
        self.assertEqual(io.write_table(outfile = self.path,
                                rows = rows,
                                header = header,
                                delimiter = delimiter,
                                verbose = verbose),
                         3)
        if not os.path.isfile(self.path):
            self.fail("File {} was not create".format(self.path))
        else:
            print("\tFile {} created".format(self.path))
        
class TestClean_dirs(unittest.TestCase):
    """Test clean dirs"""
    
    def setUp(self):
        self.path = 'temp_test_clean_dirs'
        if os.path.isdir(self.path):
            raise OSError("Directory {} existes. Cannot overwrite".format(self.path))
        os.mkdir(self.path)
        os.mkdir(self.path + "/a")
        os.mkdir(self.path + "/b")
    
    def tearDown(self):
        os.rmdir(self.path)
    
    def test_clean_dirs_list(self):
        self.assertEqual(io.clean_dirs(['a','b'], self.path),
                         ['temp_test_clean_dirs/a','temp_test_clean_dirs/b'])
    
    def test_clean_dirs_single(self):
        self.assertEqual(io.clean_dirs('a', self.path),
                         ['temp_test_clean_dirs/a'])
        self.assertEqual(io.clean_dirs('b',self.path),
                         ['temp_test_clean_dirs/b'])
     
    def test_clean_dirs_file(self):
        io.write_table(self.path + "/dirs", ['a','b'])
        self.assertEqual(io.clean_dirs(self.path + "/dirs",
                                       self.path),
                         ['temp_test_clean_dirs/a','temp_test_clean_dirs/b'])
        os.unlink(self.path + "/dirs")
        
        
if __name__  == '__main__':
    unittest.main()
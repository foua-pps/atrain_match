'''
Created on Oct 18, 2010

@author: a001696
'''

class OrrbStats():
    """Abstract class for accumulating statistics from atrain_match. (What does 
    orrb stand for?)"""
    
    def __init__(self, results_files=None):
        """Create an OrrbStats object with results from *results_files*."""
        self.results_files = results_files
        if self.results_files is not None:
            self.do_stats()
    
    def do_stats(self):
        """
        Calculate all statistics and put results in instance attributes (?).
        
        """
        if not self.results_files:
            raise RuntimeError("No results files")
    
    def printout(self):
        """Generate nice printout of the results."""
        raise NotImplementedError("The printout method should be implemented in"
                                  " a subclass of OrrbStat.")
    
    def write(self, filename):
        """Write printout to a file."""
        f = open(filename, 'w')
        for l in self.printout():
            f.write(l + '\n')
        f.close()
    


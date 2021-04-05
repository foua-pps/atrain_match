""" Module containing functions to calculate various scores """
import numpy as np

class ScoreUtils:
    def __init__(self, a, b, c, d, num=None):

        """ 
        Initialize ScoreUtils class. Takes contigency table elements 
        as arguments (see below).

                                  Observed
                 |   yes           |         no             |   Total
            -----+-----------------+------------------------+---------------
        FC   yes |      hits (a)   | false alarms (b)       | forecast yes
              no |     misses (c)  | correct negatives (d)  | forecast no
            -----+-----------------+------------------------+---------------
             Tot |  observed yes   | observed no            |    total (n)

        Argument num is an optional argument to insert a different total sum 
        for statistics using e.g. n-1 as denominator. If not provided, num 
        is calculated as a + b + c + d.

        """
        self.a = a
        self.b = b
        self.c = c
        self.d = d
        if num is not None:
            self.n = num
        else:
            self.n = a + b + c + d

    def hitrate(self):
        """
        Calculates the hitrate (fraction correct).

        Range: 0 to 1
        Perfect score: 1
        """
        return np.divide(self.a + self.d, self.n)

    def pod_0(self):
        """
        Calculates probability of detecting 0 (clr/wat).

        Range: 0 to 1
        Perfect score: 1
        """
        return np.divide(self.d, self.b + self.d)

    def pod_1(self):
        """
        Calculates probability of detecting 1 (cld/ice).

        Range: 0 to 1
        Perfect score: 1
        """
        return np.divide(self.a, self.a + self.c)

    def far_0(self):
        """
        Calculates the 0 (clr/wat) false alarm rate.

        Range: 0 to 1
        Perfect score: 0
        """
        return np.divide(self.c, self.c + self.d)

    def far_1(self):
        """
        Calculates the 1 (cld/ice) false alarm rate.

        Range: 0 to 1
        Perfect score: 0
        """
        return np.divide(self.b, self.a + self.b)

    def pofd_0(self):
        """
        Calculates the probability of false detection for 0 (clr/wat).

        Range: 0 to 1
        Perfect score: 0
        """
        return np.divide(self.c, self.a + self.c)

    def pofd_1(self):
        """
        Calculates the probability of false detection for 1 (cld/ice).

        Range: 0 to 1
        Perfect score: 0
        """
        return np.divide(self.b, self.b + self.d)

    def heidke(self):
        """
        Calculates the Heidke skill score.

        Range: -1 to 1, 0 indicates no skill
        Perfect score: 1
        """
        num = 2 * ((self.a * self.d) - (self.b * self.c))
        denom = ((self.a + self.c) * (self.c + self.d) +\
                 (self.a + self.b) * (self.b + self.d))
        return np.divide(num, denom)

    def kuiper(self):
        """
        Calculates the Hansen Kuiper skill score.

        Range: -1 to 1, 0 indicates no skill
        Perfect score: 1
        """
        num = (self.a * self.d - self.b * self.c)
        denom = (self.a + self.c) * (self.b + self.d)
        return np.divide(num, denom)

    def bias(self):
        """ Calculates the bias (mean error). """
        return np.divide(self.b - self.c, self.n)

    def mean(self):
        """ Calculates the mean. """       
        return np.divide(self.a + self.c, self.n)

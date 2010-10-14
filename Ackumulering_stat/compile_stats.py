'''
Created on Oct 13, 2010

@author: a001696
'''
import orrb_CFC_stat, orrb_CFC_stat_emissfilt, orrb_CFC_stat_surfaces
import orrb_CTH_stat, orrb_CTH_stat_emissfilt, orrb_CTH_stat_surfaces
import orrb_CTY_stat


def compile_stats():
    cfc_results = orrb_CFC_stat.do_stats()
    orrb_CFC_stat_emissfilt.do_stats()
    orrb_CFC_stat_surfaces.do_stats()
    
    orrb_CTH_stat.do_stats()
    orrb_CTH_stat_emissfilt.do_stats()
    orrb_CTH_stat_surfaces.do_stats()
    
    orrb_CTY_stat.do_stats(cfc_results)


if __name__ == '__main__':
    compile_stats()


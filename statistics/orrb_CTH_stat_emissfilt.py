# Program orrb_CTH_stat_emissfilt

# This program calculates basic statistics for the cloud top (CTH) product for
# each month

from orrb_CTH_stat import CloudTopStats as CloudTopFilteredStats


if __name__ == "__main__":
    stats = CloudTopFilteredStats()
    stats.old_interface(modes=['EMISSFILT'], output_file_desc="cth_results_emissfilt_summary")
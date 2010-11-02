# Program orrb_CTH_stat

# This program calculates basic statistics for the cloud top (CTH) product for
# each month

from orrb_CTH_stat import CloudTopStats as CloudTopSurfacesStats


if __name__ == "__main__":
    import setup
    stats = CloudTopSurfacesStats()
    stats.old_interface(modes=setup.SURFACES, output_file_desc="cth_results_surfaces_summary")
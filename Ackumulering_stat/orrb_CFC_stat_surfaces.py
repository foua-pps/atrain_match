# Program orrb_CFC_stat_surfaces

# This program calculates basic statistics for the cloud amount (CFC) product for
# each month

from orrb_CFC_stat import CloudFractionStats

# -----------------------------------------------------
class CloudFractionSurfacesStats(CloudFractionStats):
    
    def printout(self):
        lines = []
        try:
            lines.append("Month is:  %s" % self.month)
        except KeyError:
            pass
        lines.append("Total number of matched scenes is: %s" % self.scenes)
        lines.append("")
        lines.append("Total number of CALIOP matched FOVs: %d" % self.samples_cal)
        lines.append("Mean CFC CALIOP: %f" % self.mean_CFC_cal)
        lines.append("Mean error: %f" % self.bias_cal_perc)
        lines.append("RMS error: %f" % self.rms_cal)
        lines.append("Mean error MODIS: %f" % self.bias_modis_perc)
        lines.append("RMS error MODIS: %f" % self.rms_modis)
        lines.append("")
        
        for l in lines:
            print(l)
        
        return lines


if __name__ == "__main__":
    import setup
    stats = CloudFractionSurfacesStats()
    stats.old_interface(modes=setup.SURFACES, output_file_desc="cfc_results_surfaces_summary")
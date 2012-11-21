#"""
#The CloudFractionFilteredStats class calculates basic statistics for the
#cloud amount (CFC) product for each month after filtering of thin high clouds
#"""

import pdb
from orrb_CFC_stat import CloudFractionStats
class CloudFractionFilteredStats(CloudFractionStats):
    
    def printout(self):
        try:
            lines = []
            lines.append("Total number of matched scenes is: %s" % self.scenes)
            if hasattr(self, 'samples_csa'):
                lines.append("Total number of Cloudsat matched FOVs: %d" % self.samples_csa)
                lines.append("Mean CFC Cloudsat: %f" % self.mean_CFC_csa)
                lines.append("Mean error: %f" % self.bias_csa_perc)
                lines.append("RMS error: %f" % self.rms_csa)
            else:
                lines.append("No Cloudsat matched FOVs")


            lines.append("Total number of CALIOP matched FOVs: %d" % self.samples_cal)
            lines.append("Mean CFC CALIOP: %f" % self.mean_CFC_cal)
            lines.append("Mean error: %f" % self.bias_cal_perc)
            lines.append("RMS error: %f" % self.rms_cal)
            lines.append("POD cloudy: %f" % self.pod_cloudy_cal)
            lines.append("POD clear: %f" % self.pod_clear_cal)
            lines.append("FAR cloudy: %f" % self.far_cloudy_cal)
            lines.append("FAR clear: %f" % self.far_clear_cal)
            lines.append("Kuipers: %f" % self.kuipers)
            lines.append("Hitrate: %f" % self.hitrate)
            if  hasattr(self, 'bias_modis_perc'):
                lines.append("Mean error MODIS: %f" % self.bias_modis_perc)
                lines.append("RMS error MODIS: %f" % self.rms_modis)
            else:
                lines.append("No MODIS matchups")
            lines.append("")
                
        except AttributeError:
            import traceback
            traceback.print_exc()
            lines = []
            lines.append("Total number of matched scenes is: %s" % self.scenes)
            lines.append("Total number of CALIOP matched FOVs: %d" % self.samples_cal)
            lines.append("Mean CFC CALIOP: %f" % self.mean_CFC_cal)
            lines.append("Mean error: %f" % self.bias_cal_perc)
            lines.append("RMS error: %f" % self.rms_cal)
            lines.append("POD cloudy: %f" % self.pod_cloudy_cal)
            lines.append("POD clear: %f" % self.pod_clear_cal)
            lines.append("FAR cloudy: %f" % self.far_cloudy_cal)
            lines.append("FAR clear: %f" % self.far_clear_cal)
            lines.append("Kuipers: %f" % self.kuipers)
            lines.append("Hitrate: %f" % self.hitrate)
            lines.append("")


        return lines


if __name__ == "__main__":
    pdb.set_trace()
    stats = CloudFractionFilteredStats()
    stats.old_interface(modes=['EMISSFILT'], output_file_desc="cfc_results_emissfilt_summary")

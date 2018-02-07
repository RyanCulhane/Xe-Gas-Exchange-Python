# This file defines the thresholds of the binning method
mean_vent = 0.5089
std_vent = 0.1899
thre_vent = [mean_vent-2*std_vent, mean_vent-std_vent, mean_vent,
             mean_vent+std_vent, mean_vent+2*std_vent]

mean_bar = 0.4859
std_bar = 0.1484
thre_bar = [mean_bar-2*std_bar, mean_bar-std_bar, mean_bar,mean_bar+std_bar
            mean_bar+2*std_bar, mean_bar+3*std_bar, mean_bar+4*std_bar]

mean_rbc = 0.2592
std_rbc = 0.1038
thre_rbc = [mean_rbc-2*std_rbc, mean_rbc-std_rbc, mean_rbc,
             mean_rbc+std_rbc, mean_rbc+2*std_rbc]

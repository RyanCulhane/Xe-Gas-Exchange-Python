import numpy as np

# This file defines the thresholds of the binning method
mean_vent = 0.5089
std_vent = 0.1899
thre_vent = [mean_vent-2*std_vent, mean_vent-std_vent, mean_vent,
             mean_vent+std_vent, mean_vent+2*std_vent]

mean_bar = 0.4859
std_bar = 0.1484
thre_bar = [mean_bar-2*std_bar, mean_bar-std_bar, mean_bar, mean_bar+std_bar,
            mean_bar+2*std_bar, mean_bar+3*std_bar, mean_bar+4*std_bar]

mean_rbc = 0.2592
std_rbc = 0.1038
thre_rbc = [mean_rbc-2*std_rbc, mean_rbc-std_rbc, mean_rbc,
            mean_rbc+std_rbc, mean_rbc+2*std_rbc]

# bar2gas [brown color from colorbrewer]
barmap= np.zeros((9,3))
barmap[0,:] = [0, 0, 0] # background
barmap[1,:] = [1, 0, 0] # lowest value,red
barmap[2,:] = [1, 0.7143, 0]  # yellow 2
barmap[3,:] = [0.4, 0.7, 0.4] # yellow - GREEN 3
barmap[4,:] = [0, 1, 0] # GREEN 4
barmap[5,:] = [184.0/255.0, 226.0/255.0, 145.0/255.0]
barmap[6,:] = [243.0/255.0, 205.0/255.0, 213.0/255.0]
barmap[7,:] = [225.0/255.0, 129.0/255.0, 162.0/255.0]
barmap[8,:] = [197.0/255.0, 27.0/255.0, 125.0/255.0] # highest value

long_index2color = {
        1: [0, 0, 0],
        2: [1, 0, 0],
        3: [1, 0.7143, 0],
        4: [0.4, 0.7, 0.4],
        5: [0, 1, 0],
        6: [184.0/255.0, 226.0/255.0, 145.0/255.0],
        7: [243.0/255.0, 205.0/255.0, 213.0/255.0],
        8: [225.0/255.0, 129.0/255.0, 162.0/255.0],
        9: [197.0/255.0, 27.0/255.0, 125.0/255.0]
    }

 # ventilation and RBC
mmap= np.zeros((7,3))
mmap[0,:] =[0, 0, 0] # background
mmap[1,:] = [1, 0, 0] # defect 1
mmap[2,:] = [1, 0.7143, 0]# yellow 2
mmap[3,:] = [0.4, 0.7, 0.4]# yellow - GREEN 3
mmap[4,:] = [0, 1, 0]#  GREEN 4
mmap[5,:] = [0, 0.57, 0.71] # green 5
mmap[6,:] = [0, 0, 1] # high-intense 6

short_index2color = {
        1: [0, 0, 0],
        2: [1, 0, 0],
        3: [1, 0.7143, 0],
        4: [0.4, 0.7, 0.4],
        5: [0, 1, 0],
        6: [0, 0.57, 0.71],
        7: [0, 0, 1]
    }

## parameters for histograms
rbchistogram = {
        'color': (0.8941,0.1020,0.1098),
        'x_lim': 1.2,
        'y_lim': 0.1,
        'num_bins': 50,
        'refer_fit': (0.06106, 0.2604, 0.1481)
}

barhistogram = {
        'color': (0.4,0.7608,0.6471),
        'x_lim': 2.5,
        'y_lim': 0.18,
        'num_bins': 70,
        'refer_fit': (0.07006, 0.4632, 0.1953)
}

venhistogram = {
        'color': (0.4196,0.6824,0.8392),
        'x_lim': 1.0,
        'y_lim': 0.07,
        'num_bins': 50,
        'refer_fit': (0.04074, 0.494, 0.278)
}

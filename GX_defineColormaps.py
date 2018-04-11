import numpy as np

# This file defines the thresholds of the binning method
mean_vent = 0.5098
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
venhistogram = {
        'color': (0.4196,0.6824,0.8392),
        'x_lim': 1.0,
        'y_lim': 0.07,
        'num_bins': 50,
        'refer_fit': (0.04074, 0.494, 0.278),
        'xticks': [0.0, 0.5, 1.0],
        'yticks': [0.02, 0.04, 0.06],
}

barhistogram = {
        'color': (0.4,0.7608,0.6471),
        'x_lim': 2.5,
        'y_lim': 0.18,
        'num_bins': 70,
        'refer_fit': (0.07006, 0.4632, 0.1953),
        'xticks': [0.0, 1.0, 2.0, 2.5],
        'yticks': [0.05, 0.10, 0.15],
}

rbchistogram = {
        'color': (247.0/255,96.0/255,111.0/255),
        'x_lim': 1.2,
        'y_lim': 0.1,
        'num_bins': 50,
        'refer_fit': (0.06106, 0.2604, 0.1481),
        'xticks': [0.0, 0.5, 1.0, 1.2],
        'yticks': [0.05, 0.10],
}
## average and std for the reference cohort of the statics
# histogram_ticks = {
#         'ven_xticks': [0.0, 0.5, 1.0],
#         'ven_yticks': [0.02, 0.04, 0.06],
#         'bar_xticks': [0.0, 1.0, 2.0, 2.5],
#         'bar_yticks': [0.05, 0.10, 0.15],
#         'rbc_xticks': [0.0, 0.5, 1.0, 1.2],
#         'rbc_yticks': [0.05, 0.10],
# }

referece_stats = {
        'r_ven_defect_ave': '1',
        'r_ven_defect_std':'<1',
        'r_ven_low_ave':'15',
        'r_ven_low_std':'6',
        'r_ven_high_ave':'17',
        'r_ven_high_std':'10',

        'r_bar_defect_ave':'<1',
        'r_bar_defect_std':'<1',
        'r_bar_low_ave':'9',
        'r_bar_low_std':'8',
        'r_bar_high_ave':'1',
        'r_bar_high_std':'4',

        'r_rbc_defect_ave':'5',
        'r_rbc_defect_std':'3',
        'r_rbc_low_ave':'21',
        'r_rbc_low_std':'6',
        'r_rbc_high_ave':'7',
        'r_rbc_high_std':'5',

        'r_RBC2barrier_ave':'0.53',
        'r_RBC2barrier_std':'0.18'
}
################################################################################
## parameters for rat mapping
mean_vent_rat = 0.6026
std_vent_rat = 0.1895
thre_vent_rat = [mean_vent_rat-2*std_vent_rat, mean_vent_rat-std_vent_rat, mean_vent_rat,
             mean_vent_rat+std_vent_rat, mean_vent_rat+1.5*std_vent_rat]

mean_bar_rat =  33.3981
std_bar_rat = 24.8021
thre_bar_rat = [mean_bar_rat-std_bar_rat, mean_bar_rat-0.5*std_bar_rat, mean_bar_rat, mean_bar_rat+0.5*std_bar_rat,
            mean_bar_rat+1.0*std_bar_rat, mean_bar_rat+1.5*std_bar_rat, mean_bar_rat+2.0*std_bar_rat]

mean_rbc_rat = 14.7685
std_rbc_rat = 14.2638
thre_rbc_rat = [mean_rbc_rat-2.0/3.0*std_rbc_rat, mean_rbc_rat-1.0/3.0*std_rbc_rat, mean_rbc_rat,
            mean_rbc_rat+std_rbc_rat, mean_rbc_rat+2.0*std_rbc_rat]

rbchistogram_rat = {
        'color': (247.0/255,96.0/255,111.0/255),
        'x_lim': 75.0,
        'y_lim': 0.1,
        'num_bins': 100,
        'refer_fit': (0.04, mean_rbc_rat, std_rbc_rat),
        'xticks': [0, 25, 50, 75],
        'yticks': [0.05, 0.10],
}

barhistogram_rat = {
        'color': (0.4,0.7608,0.6471),
        'x_lim': 150.0,
        'y_lim': 0.18,
        'num_bins': 100,
        'refer_fit': (0.05, mean_bar_rat, std_rbc_rat),
        'xticks': [0, 50, 100,150],
        'yticks': [0.05, 0.10, 0.15],
}

venhistogram_rat = {
        'color': (0.4196,0.6824,0.8392),
        'x_lim': 1.0,
        'y_lim': 0.07,
        'num_bins': 100,
        'refer_fit': (0.03, mean_vent_rat, std_vent_rat),
        'xticks': [0.0, 0.5, 1.0],
        'yticks': [0.02, 0.04, 0.06],
}

referece_stats_rat = {
        'r_ven_defect_ave': '1',
        'r_ven_defect_std':'<1',
        'r_ven_low_ave':'15',
        'r_ven_low_std':'6',
        'r_ven_high_ave':'17',
        'r_ven_high_std':'10',

        'r_bar_defect_ave':'<1',
        'r_bar_defect_std':'<1',
        'r_bar_low_ave':'9',
        'r_bar_low_std':'8',
        'r_bar_high_ave':'1',
        'r_bar_high_std':'4',

        'r_rbc_defect_ave':'5',
        'r_rbc_defect_std':'3',
        'r_rbc_low_ave':'21',
        'r_rbc_low_std':'6',
        'r_rbc_high_ave':'7',
        'r_rbc_high_std':'5',

        'r_RBC2barrier_ave':'0.53',
        'r_RBC2barrier_std':'0.18'
}

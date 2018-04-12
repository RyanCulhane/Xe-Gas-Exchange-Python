
from GX_Twix_parser import readTwix
from scipy.stats import norm
import scipy.sparse as sps
import numpy as np
import os, pdb
from GX_Recon_utils import read_twix_hdr
from GX_Spec_classmap import NMR_TimeFit

def print_fit(fit_box,Subject_ID, key):
    if(key==1):
        print "**********//Fitting Results -- Calibration//******************"
    else:
        print "*********//Fitting Results -- Dixon Spectrum//****************"

    print "{:<8}{:<10}{:<10}{:<10}".format(Subject_ID, 'RBC','Barrier','Gas')

    for k in ['Freq','FwhmL','FwhmG','Phase']:
        v = fit_box[k]
        print "{:<8}{:<10}{:<10}{:<10}".format(k, round(v[0],2),round(v[1],2),round(v[2],2))

    v = fit_box['Area']
    print("{:<8}{:<10.2E}{:<10.2E}{:<10.2E}".format('Area',v[0],v[1],v[2]))
    print("**************************************************************")
    print("******RBC:barrier = "+str(fit_box['Area'][0]/fit_box['Area'][1]))
    print "**************************************************************"

def spect_fit(twix_cali_file, twix_dixon_file, Subject_ID):
    ## spectroscopic fitting on calibration file and dixon bonus spectrum
    ## ************************************************part 1, fit on calibration

    scans, evps = readTwix(twix_cali_file)

    data = np.asarray([x.data for x in scans])
    meas_dict,_ = read_twix_hdr(evps[2][1])

    nFids = np.shape(data)[0]
    nPts = np.shape(data)[1]
    # 200 + 1 + 20
    nSkip = 100 # skip to reach steady state
    nGas = 1 #number of dedicated gas spectra for frequency reference
    nCal = 20 #number of flipangle calibration frames following the dissolved hit

    nDis = nFids - nCal - nGas
    data_dis = data[nSkip:nDis,:]
    data_dis_ave = np.average(data_dis,axis=0)
    data_gas = data[nDis,:]

    dwell_time = float(meas_dict['alDwellTime'].split()[0])*1e-9 # unit in second
    t = np.array(range(0,nPts))*dwell_time

    ## initial fit from the calibration to determine frequency and fwhm of gas and dissolved
    gasfit = NMR_TimeFit(time_signal=data_gas, t=t, area= 1e-4, freq=-84,
                         fwhm=30, phase=0, line_boardening=0,zeropad_size=10000,method='lorenzian')
    gasfit.fit_time_signal()

    disfit = NMR_TimeFit(time_signal=data_dis_ave, t=t, area=[1,1,1],freq=[0,-700,-7400],
                         fwhmL=[250,200,30],fwhmG=[0,200,0],phase=[0,0,0],line_boardening=0,zeropad_size=np.size(t),method='voigt')

    lb = np.stack(([-np.inf,-np.inf,-np.inf],[-100,-900,-np.inf],[-np.inf,-np.inf,-np.inf],[-np.inf,-np.inf,-np.inf],[-np.inf,-np.inf,-np.inf])).flatten()
    ub = np.stack(([+np.inf,+np.inf,+np.inf],[+100,-500,+np.inf],[+np.inf,+np.inf,+np.inf],[+np.inf,+np.inf,+np.inf],[+np.inf,+np.inf,+np.inf])).flatten()
    bounds = (lb,ub)



    disfit.fit_time_signal(bounds)
    # pdb.set_trace()

    ## ************************************************part 2, fit on dixon bonus
    scans, evps = readTwix(twix_dixon_file)
    data_dixon = np.asarray(scans[-1].data)

    meas_dict,_ = read_twix_hdr(evps[2][1])
    nPts = np.size(data_dixon)
    dwell_time = float(meas_dict['alDwellTime'].split()[0])*1e-9 # unit in second
    t = np.array(range(0,nPts))*dwell_time

    dixonfit = NMR_TimeFit(time_signal=data_dixon, t=t, area=[1,1,1],freq=[0,-700,-7400],
                         fwhmL=[250,200,30],fwhmG=[0,200,0],phase=[0,0,0],line_boardening=0,zeropad_size=np.size(t),method='voigt')

    dixonfit.fit_time_signal()

    # fit again with the freq and fwhm from the calibration fitting
    area = dixonfit.area
    freq = disfit.freq - gasfit.freq
    fwhmL = disfit.fwhmL
    fwhmG = disfit.fwhmG
    phase = dixonfit.phase

    fwhmL[fwhmL<0] = 0
    # pdb.set_trace()
    assert np.prod(fwhmL>0), "There is a negative value in Dissolved fitting fwhmL"
    # assert np.prod(fwhmG>0), "There is a negative value in Dissolved fitting fwhmG"

    # set up bounds to constrain frequency and fwhm change
    lb = np.stack((0.33*area,freq-1.0,0.9*fwhmL,0.9*fwhmG-1.0,[-np.inf,-np.inf,-np.inf])).flatten()
    ub = np.stack((3.00*area,freq+1.0,1.1*fwhmL,1.1*fwhmG+1.0,[+np.inf,+np.inf,+np.inf])).flatten()

    # lb = np.stack((0.33*area,[-np.inf,-np.inf,-np.inf],0.9*fwhmL,0.9*fwhmG-1.0,[-np.inf,-np.inf,-np.inf])).flatten()
    # ub = np.stack((3.00*area,[+np.inf,+np.inf,+np.inf],1.1*fwhmL,1.1*fwhmG+1.0,[+np.inf,+np.inf,+np.inf])).flatten()
    bounds = (lb,ub)

    dixonfit = NMR_TimeFit(time_signal=data_dixon, t=t, area=area,freq=freq,
                         fwhmL=fwhmL,fwhmG=fwhmG,phase=phase,line_boardening=0,zeropad_size=np.size(t),method='voigt')

    dixonfit.fit_time_signal(bounds)

    fit_box_cali = {
        'Area': disfit.area,
        'Freq': disfit.freq,
        'FwhmL': disfit.fwhmL,
        'FwhmG': disfit.fwhmG,
        'Phase': disfit.phase,
    }
    print_fit(fit_box = fit_box_cali, Subject_ID = Subject_ID, key = 1)

    fit_box = {
        'Area': dixonfit.area,
        'Freq': dixonfit.freq,
        'FwhmL': dixonfit.fwhmL,
        'FwhmG': dixonfit.fwhmG,
        'Phase': dixonfit.phase,
    }
    print_fit(fit_box = fit_box, Subject_ID = Subject_ID, key =2)

    # pdb.set_trace()
    RBC2barrier =  dixonfit.area[0]/dixonfit.area[1]

    return RBC2barrier, fit_box

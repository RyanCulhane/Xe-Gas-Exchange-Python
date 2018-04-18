# -*- coding: utf-8 -*-
from GX_Twix_parser import readTwix
import numpy as np
import pdb, os, sys
from GX_Recon_utils import read_twix_hdr
from GX_Spec_classmap import NMR_TimeFit
from bounded_lsq.least_squares import least_squares

# twix_cali_file = 'meas_005002_cali.dat'
def spect_calibration(twix_cali_file,result_file_path):
    ## ************************************************part 1, fit on calibration
    scans, evps = readTwix(twix_cali_file)

    data = np.asarray([x.data for x in scans])
    meas_dict,_ = read_twix_hdr(evps[2][1])
    dicom_dict,_ = read_twix_hdr(evps[1][1])

    # fetch useful values
    nFids = np.shape(data)[0]
    nPts = np.shape(data)[1]
    dwell_time = float(meas_dict['alDwellTime'].split()[0])*1e-9 # unit in second
    t = np.array(range(0,nPts))*dwell_time
    te = float(meas_dict['alTE'].split()[0])
    freq = int(dicom_dict['lFrequency'])

    # 200 + 1 + 20
    nSkip = 100 # skip to reach steady state
    nGas = 1 #number of dedicated gas spectra for frequency reference
    nCal = 20 #number of flipangle calibration frames following the dissolved hit
    nDis = nFids - nGas - nCal

    # some parameters that may be useful
    freq_std = 34091550
    freq_tol = 200
    deltaPhase1_tol = 90
    TrueRefScale_tol = 1.17
    SNR_tol = 25
    flip_angle_target = 20
    rbc_bar_adjust = 80.0

    # split gas and dissolved data
    nDis = nFids - nCal - nGas
    data_dis = data[nSkip:nDis,:]
    data_dis_ave = np.average(data_dis,axis=0)
    data_gas = data[nDis,:]

    ## initial fit from the calibration to determine frequency and fwhm of gas and dissolved
    gasfit = NMR_TimeFit(time_signal=data_gas, t=t, area= 1e-4, freq=-84,
                         fwhm=30, phase=0, line_boardening=0,zeropad_size=10000,method='lorenzian')
    gasfit.fit_time_signal()

    disfit = NMR_TimeFit(time_signal=data_dis_ave, t=t, area=[1,1,1],freq=[0,-700,-7400],
                         fwhmL=[250,200,30],fwhmG=[0,200,0],phase=[0,0,0],line_boardening=0,zeropad_size=np.size(t),method='voigt')
    disfit.fit_time_signal()

    # 1 . report frequency
    freq_target = freq + gasfit.freq
    freq_target = int(np.asscalar(np.round(freq_target)))

    print('\r\nFrequency_target = {:8.0f} Hz'.format(freq_target))
    stream = 'Frequency_target = {:8.0f} Hz *****\r\n'.format(freq_target)

    if((freq_target - freq_std) > freq_tol):
        print('***Warning! Frequency adjust exceeds tolerances; Check system')
        stream = stream + '***Warning! Frequency adjust exceeds tolerances; Check system\r\n'

    # 2. report TE90
    deltaPhase = disfit.phase[1] - disfit.phase[0]
    deltaPhase = np.mod(abs(deltaPhase),180)
    deltaFreq = abs(disfit.freq[1] - disfit.freq[0])
    deltaTE90 = (90 - deltaPhase)/(360*deltaFreq)
    TE90 = te + deltaTE90*1e6 # in usec

    print("TE90 = {:3.2f} ms".format(TE90/1000))
    stream = stream + "TE90 = {:3.2f} ms *****\r\n".format(TE90/1000)

    if abs(deltaPhase)> 90 :
        print('***WARNING! Phi_cal = {:3.0f}{}; Use min TE!'.format(deltaPhase,u'\u00b0'.encode('utf8')))
        stream = stream + '***WARNING! Phi_cal = {:3.0f}{}; Use min TE!\r\n'.format(deltaPhase,u'\u00b0'.encode('utf8'))

    # 3. report reference voltage (Flip angle)
    calData = data[(nDis+1):(nDis+nCal+1),:]
    flipCalAmps = np.amax(abs(calData),axis=1)

    guess = [np.max(flipCalAmps), 20.0*np.pi/180]

    x_data = np.array(range(1,len(flipCalAmps)+1))
    y_data = flipCalAmps

    # curve fitting using trust region reflection algorithm
    def calc_fitting_residual(coef):
        y_fit = coef[0]*np.cos(coef[1])**(x_data - 1)
        residual = (y_fit - y_data).flatten()
        return residual

    max_nfev = 13000
    bounds = ([-np.inf, -np.inf],[np.inf, np.inf])

    fit_result = least_squares(fun=calc_fitting_residual, x0=guess, jac='2-point', bounds=bounds,
                               method='dogbox', max_nfev=max_nfev)

    flip_angle = abs(fit_result['x'][1]*180/np.pi)
    v_ref_scale_factor = flip_angle_target/flip_angle

    print('True Ref scale factor = {:3.3f}'.format(v_ref_scale_factor))
    stream = stream + 'True Ref scale factor = {:3.3f} *****\r\n'.format(v_ref_scale_factor)

    print('For 600V calibration, True_Ref = {:3.0f} V'.format(600*v_ref_scale_factor))
    stream = stream +'For 600V calibration, True_Ref = {:3.0f} V *****\r\n'.format(600*v_ref_scale_factor)

    if(v_ref_scale_factor > TrueRefScale_tol):
        print('***Warning! Excessive calibration scale factor; check system')
        stream = stream + '***Warning! Excessive calibration scale factor; check system\r\n'

    # 4. calculate and report SNR for dissolved peaks
    time_fit = disfit.calc_time_sig(disfit.t)
    time_res = time_fit - disfit.time_signal
    n25pct = int(round(len(time_res)/4.0))
    std25 = np.std(time_res[-n25pct:])
    SNR_dis = disfit.area/std25
    print('\r\nSNR for dissolved peaks = {:3.1f}, {:3.1f}, {:3.1f}'.format(SNR_dis[0], SNR_dis[1], SNR_dis[2]))
    stream = stream + '\r\nSNR for dissolved peaks = {:3.1f}, {:3.1f}, {:3.1f}\r\n'.format(SNR_dis[0], SNR_dis[1], SNR_dis[2])

    # 5. calculate and report SNR for gas
    time_fit = gasfit.calc_time_sig(gasfit.t)
    time_res = time_fit - gasfit.time_signal
    n25pct = int(round(len(time_res)/4.0))
    std25 = np.std(time_res[-n25pct:])
    SNR_gas=gasfit.area/std25
    print('SNR for gas peak = {:3.1f}'.format(SNR_gas[0]))
    stream = stream + 'SNR for gas peak = {:3.1f}\r\n'.format(SNR_gas[0])

    if SNR_gas < SNR_tol:
        print('***WARNING! Gas FID SNR below minimums; Check Coil Plug')
        stream = stream + '***WARNING! Gas FID SNR below minimums; Check Coil Plug\r\n'

    # 6. Quantify ammount of off resonance excitation
    gas_dis_ratio = disfit.area[2]/sum(disfit.area[:2])
    print('\r\ngas_dis_ratio = {:3.3f}'.format(gas_dis_ratio))
    stream = stream + '\r\ngas_dis_ratio = {:3.3f}\r\n'.format(gas_dis_ratio)

    # 7. Quantify RBC:Barrier ratio
    rbc_bar_ratio = disfit.area[0]/disfit.area[1]
    print('rbc_bar_ratio = {:3.3f}'.format(rbc_bar_ratio))
    print('RbcBar_{:2.0f} = {:3.3f}'.format(rbc_bar_adjust,rbc_bar_ratio*rbc_bar_adjust/100.0))

    stream = stream + 'rbc_bar_ratio = {:3.3f}\r\n'.format(rbc_bar_ratio)
    stream = stream + 'RbcBar_{:2.0f} = {:3.3f}\r\n'.format(rbc_bar_adjust,rbc_bar_ratio*rbc_bar_adjust/100.0)

    with open(result_file_path,'w') as cali_doc:
        cali_doc.write(stream)

# main function here
if __name__ == "__main__":

    if (len(sys.argv) == 2):
        twix_cali_file = sys.argv[1]
        result_file_path = 'calibration.txt'

    elif(len(sys.argv) == 3):
        twix_cali_file = sys.argv[1]
        result_file_path = sys.argv[2]

    elif(len(sys.argv) == 1):
        twix_cali_file = 'meas_005002_cali.dat'
        result_file_path = 'calibration_result.txt'

    else:
        print "Usage 1: python GX_Spec_calibration.py <calibration twix file> <output file direction/name>"
        print "Usage 1: python GX_Spec_calibration.py <calibration twix file>"
        sys.exit(-1)

    spect_calibration(twix_cali_file = twix_cali_file,result_file_path=result_file_path)

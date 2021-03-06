# scrip auto check if there is a scan completed, fetch the data, and launch the gas exhcange mapping program

import os,re
import pdb
import glob
import datetime as dt

from  GX_Spec_calibration import spect_calibration, send_sms

current_day = dt.date.today()
# current_day = dt.date(2018,3,5)

# find all folders that are modified today
drive = '/media/rawdata'
cali_txt_name = 'calibration_results.txt'

# check connection
if os.listdir(drive) == []:
    print('Drive empty, probabely lost connection, check time: '+ dt.datetime.now().strftime('%Y-%m-%d %H:%M:%S'))
else:
    # acquire list of folders
    drive_subjects =  [ os.path.join(drive,folder) for folder in os.listdir(drive) if os.path.isdir(os.path.join(drive, folder))]

    for path in drive_subjects:

        st = os.stat(path)

        # folders that are created or modiftied within 3 days
        if((current_day - dt.date.fromtimestamp(st.st_mtime)) < dt.timedelta(3)):

            list_cali = glob.glob(path+'/meas*Xe_fid_cali*.dat')
            key_exist = os.path.isfile(os.path.join(path,cali_txt_name))

            if((len(list_cali)>0) & (~key_exist)):

                ## parse subject name
                subject_p = re.compile('[0-9]{3}-[0-9]{3}[A-F]?')
                subject_name = subject_p.findall(path)
                if(len(subject_name)>0):
                    subject_name = subject_name[0]
                else:
                    subject_name = ''

                ## call calibration function
                twix_cali_file = os.path.join(path,list_cali[0])
                result_file_path = os.path.join(path,cali_txt_name)

                stream = spect_calibration(twix_cali_file = twix_cali_file,result_file_path = result_file_path)
                
                #print(subject_name)
                send_sms(stream+'This is for Subejct '+subject_name+'\r\n')
                print('Calibration for done in: ' + path)

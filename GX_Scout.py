# scrip auto check if there is a scan completed, fetch the data, and launch the gas exhcange mapping program

import os
import pdb
import glob
import datetime as dt
import re
from shutil import copy2

current_day = dt.date.today()
# current_day = dt.date(2018,3,5)

# find all folders that are modified today
local = '/home/ziyiw/Patients'
drive = '/media/rawdata'

# acquire list of subjects that have been processed
subject_pattern = re.compile('[0-9]{6}')
drive_subject_p = re.compile('[0-9]{3}-[0-9]{3}[A-F]?')
key = 0

if os.listdir(drive) == []:
    print('Drive empty, probabely lost connection, check time: '+ dt.datetime.now().strftime('%Y-%m-%d %H:%M:%S'))
else:
    # check, if the driver is connected
    local_subjects = [ name for name in os.listdir(local) if os.path.isdir(os.path.join(local, name)) ]

    # this command will match for exact name
    local_subjects = [ name for name in local_subjects if re.match(subject_pattern, name)]

    # acquire list of new subject from the driver
    drive_subjects =  [ os.path.join(drive,item) for item in os.listdir(drive) if os.path.isdir(os.path.join(drive, item))]

    for path in drive_subjects:

        st = os.stat(path)

        # folders that are created or modiftied within 3 days
        if((current_day - dt.date.fromtimestamp(st.st_mtime)) < dt.timedelta(3)):

            list_cali = glob.glob(path+'/meas*Xe_fid_cali*.dat')
            list_ute = glob.glob(path+'/meas*BHUTE_Radial*.dat')
            list_dixon = glob.glob(path+'/meas*Xe_Radial_Dixon*.dat')

            # if all three files exit
            if(len(list_cali)*len(list_ute)*len(list_dixon) > 0):

                # using regular expression pattern match to find subject name
                subject_name = drive_subject_p.findall(path)

                if len(subject_name) != 1:
                    print("new potential scan found, illegal naming!: "+path)
                    continue

                subject_name = subject_name[0]
                subject_name = subject_name[0:3]+subject_name[4:]
                pdb.set_trace()

                # ensure that the subject has not been previously processed
                if (subject_name not in local_subjects):

                    print('new scan completed: '+subject_name)
                    key = 1

                    subject_path = local+'/'+subject_name

                    os.makedirs(subject_path)

                    # copy files to local
                    copy2(list_cali[0],subject_path)
                    copy2(list_dixon[0],subject_path)
                    copy2(list_ute[0],subject_path)

                    # lauch the process program
                    from GX_Process import GX_Process
                    GX_Process(data_dir=local+'/'+subject_name, Subject_ID=subject_name)
                    # finish the process program
                    print('Process finished, complete time: '+ dt.datetime.now().strftime('%Y-%m-%d %H:%M:%S'))

                # in case we ran a second Dixon
                if(len(list_dixon) > 1):

                    subject_name = subject_name +'_highBW'

                    if (subject_name not in local_subjects):

                        print('A second new scan completed: '+subject_name)
                        key = 1

                        subject_path = local+'/'+subject_name

                        os.makedirs(subject_path)

                        # copy files to local
                        copy2(list_cali[0],subject_path)
                        copy2(list_dixon[1],subject_path)
                        copy2(list_ute[-1],subject_path)

                        # lauch the process program
                        from GX_Process import GX_Process
                        GX_Process(data_dir=local+'/'+subject_name, Subject_ID=subject_name)
                        # finish the process program
                        print('Process finished, complete time: '+ dt.datetime.now().strftime('%Y-%m-%d %H:%M:%S'))


    if(key==0):
        print('No new subject found, check time: '+ dt.datetime.now().strftime('%Y-%m-%d %H:%M:%S'))

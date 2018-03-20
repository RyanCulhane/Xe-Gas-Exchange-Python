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
local = '/home/ziyi/Gas_Exchange'
drive = '/media/ziyi'

# acquire list of subjects that have been processed
subject_pattern = re.compile('[0-9]{6}')

local_subjects = [ name for name in os.listdir(local) if os.path.isdir(os.path.join(local, name)) ]

local_subjects = [ name for name in local_subjects if re.match(subject_pattern, name)]

# acquire list of new subject from the driver
drive_subjects =  [ os.path.join(drive,item) for item in os.listdir(drive) if os.path.isdir(os.path.join(drive, item))]

for path in drive_subjects:

    st = os.stat(path)

    # folders that are created or modiftied today
    if(current_day == dt.date.fromtimestamp(st.st_mtime)):

        list_cali = glob.glob(path+'/meas*Xe_fid_cali*.dat')
        list_ute = glob.glob(path+'/meas*1H_BHUTE_Radial*.dat')
        list_dixon = glob.glob(path+'/meas*Xe_Radial_Dixon*.dat')

        # if all three files exit
        if(len(list_cali)*len(list_ute)*len(list_dixon) > 0):

            subject_name = path[-7:-4] + path[-3:]

            # ensure that the subject has not been previously processed
            if (subject_name not in local_subjects):

                print('new scan completed: '+subject_name)

                os.makedirs(subject_name)

                # copy files to local
                copy2(list_cali[0],subject_name)
                copy2(list_dixon[0],subject_name)
                copy2(list_ute[0],subject_name)

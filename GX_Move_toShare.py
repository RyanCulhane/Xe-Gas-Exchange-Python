# this script move all the localy processed color stack into the share folder
# Ziyi. 03.29.2018

import os
import glob
from shutil import copyfile


target_p = '/home/ziyiw/Patients_old/'
destination_p = '/media/sharedrive/shared/_ziyi/ColorStacks_Siemens/'
des_rep = '/media/sharedrive/shared/_ziyi/Report_Siemens/'

des_ven = destination_p+'VentilationCor/'
des_bar = destination_p+'Barrier2gas/'
des_rbc = destination_p+'RBC2gas/'

subject_list = os.listdir(target_p)

for subject in subject_list:
    file_name = 'ven_Sub'+subject+'.nii'
    copyfile(target_p+subject+'/'+file_name, des_ven+file_name)

    file_name = 'Bar2gas_Sub'+subject+'.nii'
    copyfile(target_p+subject+'/'+file_name, des_bar+file_name)

    file_name = 'RBC2gas_Sub'+subject+'.nii'
    copyfile(target_p+subject+'/'+file_name, des_rbc+file_name)

    file_name = 'report_'+subject+'.pptx'
    copyfile(target_p+subject+'/'+file_name, des_rep+file_name)

    print('completed transition of Subeject '+subject)

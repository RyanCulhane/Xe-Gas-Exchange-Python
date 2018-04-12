import sys

def GX_Process_rat(data_dir, Subject_ID):

    from GX_ratmap_classmap import GXRat
    # from GX_email import send_email

    # start = time.time()
    rat = GXRat(data_dir=data_dir, Subject_ID=Subject_ID)
    rat.makeReport()
    rat.saveMat()

if __name__ == "__main__":

    # Specify where is your files (grouped by folder)
    data_dir = '/media/sharedrive/shared/team_documents/2018_Preclinical_PAH//RawData/HEALTHY/'

    # Specify what subjects you need to process
    Subject_IDs = ['PAH1','PAH2','PAH3','PAH4','PAH5','PAH6','PAH7','PAH8','PAH9','PAH10']
    # Subject_IDs = ['PAH7B','PAH10B','PAH11','PAH12']
    # Subject_IDs = ['PAH13','PAH14']

    for Subject_ID in Subject_IDs:
        file_dir = data_dir+Subject_ID
        GX_Process_rat(data_dir=file_dir, Subject_ID=Subject_ID)
        print('Finished Rat: '+Subject_ID)

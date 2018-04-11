import sys

def GX_Process_rat(data_dir, Subject_ID):

    from GX_ratmap_classmap import GXRat
    # from GX_email import send_email

    # start = time.time()
    rat = GXRat(data_dir=data_dir, Subject_ID=Subject_ID)
    rat.makeReport()
    rat.saveMat()

    # send out report via email
    # send_email(data_dir=subject.data_dir, Subject_ID=subject.Subject_ID)

if __name__ == "__main__":

    # if (len(sys.argv) == 2):
    #     # data_dir = '/media/sharedrive/shared/_ziyi/Patients/'+sys.argv[1]
    #     data_dir = '/media/sharedrive/shared/team_documents/2018_Preclinical_PAH//RawData/HEALTHY/'+sys.argv[1]
    #     Subject_ID = sys.argv[1]
    #
    # elif(len(sys.argv) == 3):
    #     data_dir = sys.argv[1]
    #     Subject_ID = sys.argv[2]
    #
    # else:
    #     print "Usage 1: python GX_classmap.py <data directory/Subject_ID>"
    #     print "Usage 2: python GX_classmap.py <data directory> <Subject_ID>"
    #     sys.exit(-1)

    # Subject_IDs = ['PAH2']
    Subject_IDs = ['PAH1','PAH2','PAH3','PAH4','PAH5','PAH6','PAH7','PAH8','PAH9','PAH10']

    for Subject_ID in Subject_IDs:
        data_dir = '/media/sharedrive/shared/team_documents/2018_Preclinical_PAH//RawData/HEALTHY/'+Subject_ID
        GX_Process_rat(data_dir, Subject_ID)
        print('Finished Rat: '+Subject_ID)

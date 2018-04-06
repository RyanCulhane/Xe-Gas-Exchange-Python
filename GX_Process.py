import sys

def GX_Process(data_dir, Subject_ID):

    from GX_Map_classmap import GXSubject

    # start = time.time()
    subject = GXSubject(data_dir=data_dir, Subject_ID=Subject_ID)
    # subject.readFromMat()
    subject.GXRecon()
    subject.GXMapping()
    subject.saveMat()
    # subject.sendEmail()

if __name__ == "__main__":

    if (len(sys.argv) == 2):
        # data_dir = '/media/sharedrive/shared/_ziyi/Patients/'+sys.argv[1]
        data_dir = '/home/ziyiw/Patients/'+sys.argv[1]
        Subject_ID = sys.argv[1]

    elif(len(sys.argv) == 3):
        data_dir = sys.argv[1]
        Subject_ID = sys.argv[2]

    else:
        print "Usage 1: python GX_classmap.py <data directory/Subject_ID>"
        print "Usage 2: python GX_classmap.py <data directory> <Subject_ID>"
        sys.exit(-1)

    GX_Process(data_dir, Subject_ID)

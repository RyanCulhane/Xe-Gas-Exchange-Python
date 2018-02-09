import numpy as np
import numpy as np
import nibabel as nib
import scipy.io as sio
import sys

from matplotlib import pyplot as plt

import pdb

class GXSubject(object):
    def __init__(self):
        print("Initiate subject 002102. Today is a beautiful day!")
        self.filename = []
        self.RBC2barrier = 0.468
        self.subjectID = []
        self.TE90 = 460

        self.gas_highreso = []
        self.gas_highSNR = []
        self.dissolved = []

        self.rbc = []
        self.barrier = []
        self.gas_binning = []
        self.rbc2gas_binning = []
        self.barrier2gas_binning = []

        self.ute = []
        self.mask = []
        self.ute_reg = []
        self.mask_reg = []
        self.mask_reg_vent = []

        self.key_box = {}
        self.stats = {}

        print("read in Xe data")
        self.readinXe()

        print("mask processing")
        self.uteSegmentation()
        # print("check alignemnt")
        # self.checkAlignment()
        self.uteRegister()

        print("Gas_highreso binning and mask_vent")
        self.gasBinning()

        print("Dixon and binning")
        self.dixonDecomp()
        self.barBinning()
        self.rbcBinning()

        print("Clnical Report")

    def checkAlignment(self):
        # check if all the images are the same alignment
        dim1 = 80
        dim2 = 50
        plt.figure()
        plt.subplot(2,4,1)
        plt.imshow(abs(self.gas_highreso[:,:,dim1]))
        plt.title("gas_highreso")
        plt.subplot(2,4,2)
        plt.imshow(self.mask[:,:,dim1])
        plt.title("mask")
        plt.subplot(2,4,3)
        plt.imshow(self.ute[:,:,dim1])
        plt.title("ute")
        plt.subplot(2,4,4)
        plt.imshow(abs(self.dissolved[:,:,dim1]))
        plt.title("dissolved")

        plt.subplot(2,4,5)
        plt.imshow(abs(self.gas_highreso[:,:,dim2]))
        plt.title("gas_highreso")
        plt.subplot(2,4,6)
        plt.imshow(self.mask[:,:,dim2])
        plt.title("mask")
        plt.subplot(2,4,7)
        plt.imshow(self.ute[:,:,dim2])
        plt.title("ute")
        plt.subplot(2,4,8)
        plt.imshow(abs(self.dissolved[:,:,dim2]))
        plt.title("dissolved")

        plt.show()

    def readinXe(self):
        ## temporal usage
        # read in dissolved and gas Xe
        fdata = 'Sub002102_data.mat'

        mat_input = sio.loadmat(fdata)
        self.gas_highreso = mat_input['gasVol_highreso']
        self.gas_highSNR = mat_input['gasVol_highSNR']
        self.dissolved = mat_input['dissolvedVol']

    def uteSegmentation(self):
        ## temporal usage
        from GX_CNNpredict import CNNpredict

        # fmask = 'BHUTE_Sub002102_FID49886_mask_grow.nii'
        # self.mask = np.array(nib.load(fmask).get_data(),dtype='bool')

        fute = 'BHUTE_Sub002102_FID49886_recon.nii'

        self.ute = np.array(nib.load(fute).get_data())

        self.mask = CNNpredict(ute = self.ute)

        pdb.set_trace()

    def uteRegister(self):
        ## temporal usage
        from GX_utils import register

        self.ute_reg, self.mask_reg = register(gas_highreso = abs(self.gas_highreso),
                                               ute          = self.ute,
                                               mask         = self.mask.astype(float))

    def gasBinning(self):
        ## Binning for gas_highreso
        from GX_defineColormaps import thre_vent
        from GX_utils import gasBinning

        self.gas_binning,self.mask_reg_vent = gasBinning(gas_highreso  = abs(self.gas_highreso),
                                                         bin_threshold = thre_vent,
                                                         mask          = self.mask_reg,
                                                         percentile    = 0.99)

    def dixonDecomp(self):
        ## Dixon decomposition to get rbc and barrier from dissolved
        from GX_utils import dixonDecomp

        self.rbc, self.barrier = dixonDecomp(gas_highSNR     = abs(self.gas_highSNR),
                                             dissolved       = self.dissolved,
                                             mask_vent       = self.mask_reg_vent,
                                             meanRbc2barrier = self.RBC2barrier)
    def barBinning(self):
        ## binning for barrier
        from GX_defineColormaps import thre_bar
        from GX_utils import disBinning

        cor_TE90 = np.exp(self.TE90/2000)/np.exp(self.TE90/50000)
        cor_flipoff = 100*np.sin(0.5*np.pi/180)/np.sin(20*np.pi/180)

        self.barrier2gas_binning = disBinning(discomp    = self.barrier,
                                           gas_highSNR   = abs(self.gas_highSNR),
                                           bin_threshold = thre_bar,
                                           mask          = self.mask_reg_vent,
                                           cor           = cor_TE90*cor_flipoff)
    def rbcBinning(self):
        ## binning for barrier
        from GX_defineColormaps import thre_rbc
        from GX_utils import disBinning

        cor_TE90 = np.exp(self.TE90/2000)/np.exp(self.TE90/50000)
        cor_flipoff = 100*np.sin(0.5*np.pi/180)/np.sin(20*np.pi/180)

        self.rbc2gas_binning = disBinning(discomp       = self.rbc,
                                          gas_highSNR   = abs(self.gas_highSNR),
                                          bin_threshold = thre_rbc,
                                          mask          = self.mask_reg_vent,
                                          cor           = cor_TE90*cor_flipoff)

if __name__ == "__main__":
    # if len(sys.argv) != 4:
    #     print "Usage: python geohash2kml.py <mode> <input file> <output file>"
    #     sys.exit(-1)
    #
    # output_type = sys.argv[1]
    # if output_type not in OUTPUT_TYPES:
    #     print "Please select one of the following modes: {}".format(OUTPUT_TYPES)
    #     sys.exit(1)

    # Collect file arguments
    # input_file = sys.argv[2]
    # output_file = sys.argv[3]

    # Create helper object
    subject = GXSubject()
    pdb.set_trace()

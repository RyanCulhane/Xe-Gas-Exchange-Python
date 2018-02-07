import numpy as np
import numpy as np
import nibabel as nib
import scipy.io as sio
import sys

from matplotlib import pyplot as plt

import pdb

class GXSubject(object):
    def __init__(self,filename,RBC2barrier,subjectID):
        print("Initiate subject "+subjectID+" Today is a beautiful day!")
        self.filename = filename
        self.RBC2barrier = RBC2barrier
        self.subjectID = subjectID
        self.TE90 = 460

        self.gas_highreso = []
        self.gas_highSNR = []
        self.dissolved = []
        self.rbc = []
        self.barrier = []
        self.gas_binning = []
        self.rbc2gas_binning = []
        self.barrier2gas_binning = []

        self.mask = []
        self.mask_reg = []
        self.mask_reg_vent = []

        self.key_box = {}
        self.stats = {}

        print("mask processing")
        self.maskInit(filename)
        self.maskRegister(filename)

        print("Gas_highreso binning and mask_vent")
        self.gasBinning()

        print("Dixon and binning")
        self.dixonDecomp()
        self.barBinning()
        self.rbcBinning()

        print("Clnical Report")

    def maskInit(self):
        ## temporal usage
        fmask = 'BHUTE_Sub002102_FID49886_mask_grow.nii'

        mask = np.array(nib.load(fmask).get_data(),dtype='bool')
        self.mask = mask

    def maskRegister(self):
        ## temporal usage
        self.mask_reg = mask_reg

    def gasBinning(self):
        ## Binning for gas_highreso
        from GX_defineColormaps import thre_vent
        from GX_utils import gasBinning

        self.gas_binning,self.mask_reg_vent = gasBinning(gas_highreso  = self.gas_highreso,
                                                         bin_threshold = thre_vent,
                                                         mask          = self.mask_reg,
                                                         percentile    = 0.99)

    def dixonDecomp(self):
        ## Dixon decomposition to get rbc and barrier from dissolved
        from GX_utils import dixonDecomp

        self.rbc, self.barrier = dixonDecomp(gas_highSNR     = self.gas_highSNR,
                                             dissolved       = self.dissolved,
                                             mask_vent       = self.mask_reg_vent,
                                             meanRbc2barrier = self.RBC2barrier)
     def barBinning(self):
         ## binning for barrier
         from GX_defineColormaps import thre_bar
         from GX_utils import disBinning

         cor_TE90 = np.exp(self.TE90/2000)/np.exp(self.TE90/50000)
         cor_flipoff = 100*np.sin(0.5*np.pi/180)/np.sin(20*np.pi/180)

         self.barrier2gas_binning = disBinning(discomp       = self.barrier,
                                               gas_highSNR   = self.gas_highSNR,
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
                                          gas_highSNR   = self.gas_highSNR,
                                          bin_threshold = thre_bar,
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

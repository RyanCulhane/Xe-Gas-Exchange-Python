import numpy as np
import numpy as np
import nibabel as nib
import scipy.io as sio
import glob, os
import sys
import time

from matplotlib import pyplot as plt

import pdb

class GXRat(object):
    def __init__(self, data_dir, Subject_ID):
        print("*********************Initiate a new rat. Today is a beautiful day!")
        self.data_dir = data_dir
        self.filename = []
        self.RBC2barrier = 0.580
        self.Subject_ID = Subject_ID
        self.TE90 = 210
        self.FOV = 5.0

        self.gas_highreso = []
        self.gas_highSNR = []
        self.gas_highreso_cor = []
        self.gas_biasfield = []
        self.dissolved = []

        self.rbc = []
        self.barrier = []
        self.ventilation = []
        self.rbc2gas = []
        self.bar2gas = []
        self.ven_binning = []
        self.rbc2gas_binning = []
        self.bar2gas_binning = []

        self.ute = []
        self.mask = []
        self.ute_reg = []
        self.mask_reg = []
        self.mask_reg_vent = []

        self.key_box = {}
        self.stats_box = {}

        print("*********************Read in 129Xe data")
        self.readinXe()

        print("*********************Mask Segmentation")
        # self.uteSegmentation()
        self.alignImages()
        # self.checkAlignment()
        self.uteRegister()

        print("*********************Gas BiasFied Correction, Binning and Ventilation Mask")
        self.gasBiasFieldCor()
        self.gasBinning()

        print("*********************Dixon and binning")
        self.dixonDecomp()
        self.barBinning()
        self.rbcBinning()

        self.generateStats()

    def makeReport(self):
        print("*********************Clnical Report")
        self.generateFigures()
        self.generateHtmlPdf()

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
        # fdata = self.data_dir+'/Sub'+self.Subject_ID+'_data.mat'

        fdata = self.data_dir+'/'+self.Subject_ID+'.mat'

        # pdb.set_trace()

        mat_input = sio.loadmat(fdata)

        # unwrape the dict
        mat_input = mat_input[self.Subject_ID]

        self.gas_highreso = mat_input['gas_highreso'][0,0]
        self.gas_highSNR = mat_input['gas_highSNR'][0,0]
        self.dissolved = mat_input['dissolved'][0,0]
        self.ute = abs(mat_input['ute'][0,0])
        self.mask = mat_input['mask'][0,0].astype(bool)
        self.airways = mat_input['airways'][0,0].astype(bool)
        self.RBC2barrier = mat_input['RBC2barrier'][0,0].flatten()
        self.TE90 = mat_input['TE90'][0,0].flatten()

    def uteSegmentation(self):
        ## temporal usage
        from GX_CNNpredict import CNNpredict

        # fmask = 'BHUTE_Sub002102_FID49886_mask_grow.nii'
        # self.mask = np.array(nib.load(fmask).get_data(),dtype='bool')
        # fute = glob.glob(self.data_dir+'/BHUTE_Sub'+self.Subject_ID+'_FID*recon.nii')[0]
        # fute = 'BHUTE_Sub002102_FID49886_recon.nii'
        # self.ute = np.array(nib.load(fute).get_data())

        # self.mask = CNNpredict(ute = self.ute)
        fmask = glob.glob(self.data_dir+'/BHUTE_Sub'+self.Subject_ID+'_FID*mask_grow.nii')[0]
        self.mask = np.array(nib.load(fmask).get_data())

    def alignImages(self):

        def alignrot(x):
            # return(np.flip(np.flip(np.flip(x,0),1),2))
            return(np.flip(np.transpose(x,(2,0,1)),2))

        self.ute = alignrot(self.ute)
        self.mask = alignrot(self.mask)
        self.gas_highreso = alignrot(self.gas_highreso)
        self.gas_highSNR = alignrot(self.gas_highSNR)
        self.dissolved = alignrot(self.dissolved)
        self.airways = alignrot(self.airways)

    def uteRegister(self):
        ## temporal usage
        from GX_Map_utils import register

        self.ute_reg, self.mask_reg = register(gas_highreso = abs(self.gas_highreso),
                                               ute          = self.ute,
                                               mask         = self.mask.astype(float))

        # take out the airways after uteRegister
        self.mask_reg = self.mask_reg & (~self.airways)

        # fmask_reg = 'BHUTE_Sub002102_FID49886_mask_grow_reg.nii'
        # self.mask_reg = np.array(nib.load(fmask_reg).get_data()).astype(bool)
        #
        # fute_reg = 'BHUTE_Sub002102_FID49886_recon_reg.nii'
        # self.ute_reg = np.array(nib.load(fute_reg).get_data())
        #
        # def alignrot(x):
        #     return(np.flip(np.flip(np.flip(x,0),1),2))
        #
        # self.ute_reg = alignrot(self.ute_reg)
        # self.mask_reg = alignrot(self.mask_reg)


    def gasBiasFieldCor(self):
        # conduct bias field correction for gas_highreso
        from GX_Map_utils import biasFieldCor

        self.gas_highreso_cor, self.gas_biasfield = biasFieldCor(image = abs(self.gas_highreso), mask = self.mask_reg)

    def gasBinning(self):
        ## Binning for gas_highreso
        from GX_defineColormaps import thre_vent_rat
        from GX_Map_utils import gasBinning

        self.ventilation, self.ven_binning, self.mask_reg_vent = gasBinning(gas_highreso  = abs(self.gas_highreso_cor),
                                                                            bin_threshold = thre_vent_rat,
                                                                            mask          = self.mask_reg,
                                                                            percentile    = 99)

    def dixonDecomp(self):
        ## Dixon decomposition to get rbc and barrier from dissolved
        from GX_Map_utils import dixonDecomp

        self.rbc, self.barrier = dixonDecomp(gas_highSNR     = self.gas_highSNR,
                                             dissolved       = self.dissolved,
                                             mask_vent       = self.mask_reg_vent,
                                             meanRbc2barrier = self.RBC2barrier)
    def barBinning(self):
        ## binning for barrier
        from GX_defineColormaps import thre_bar_rat
        from GX_Map_utils import disBinning

        # T2*_gas = 2.01 ms; T2*_rbc = 0.73 ms, T2*_bar = 0.56 ms
        # Flip angle_gas = 15, flip angle_dissolved = 20
        cor_TE90 = np.exp(self.TE90/560.0)/np.exp(self.TE90/2010.0)
        cor_flipoff = 100*np.sin(15.0*np.pi/180)/np.sin(20.0*np.pi/180)
        #
        self.bar2gas, self.bar2gas_binning = disBinning(discomp       = self.barrier,
                                                        gas_highSNR   = abs(self.gas_highSNR),
                                                        bin_threshold = thre_bar_rat,
                                                        mask          = self.mask_reg_vent,
                                                        cor           = cor_TE90*cor_flipoff)

        # from GX_defineColormaps import thre_vent
        # from GX_Map_utils import gasBinning
        # bar2gas = np.zeros(np.shape(self.barrier))
        # mask = self.mask_reg_vent
        # bar2gas[mask] = np.divide(self.barrier[mask],abs(self.gas_highSNR[mask]))
        # bar2gas[bar2gas<0] = 1e-5
        # self.bar2gas = bar2gas*cor_TE90*cor_flipoff
        #
        # _, self.bar2gas_binning,_ = gasBinning(gas_highreso  = bar2gas,
        #                                        bin_threshold = thre_vent,
        #                                        mask          = self.mask_reg_vent,
        #                                        percentile    = 99)

    def rbcBinning(self):
        ## binning for barrier
        from GX_defineColormaps import thre_rbc_rat
        from GX_Map_utils import disBinning

        # T2*_gas = 2.01 ms; T2*_rbc = 0.73 ms, T2*_bar = 0.56 ms
        # Flip angle_gas = 15, flip angle_dissolved = 20
        cor_TE90 = np.exp(self.TE90/730.0)/np.exp(self.TE90/2010.0)
        cor_flipoff = 100*np.sin(15.0*np.pi/180)/np.sin(20.0*np.pi/180)
        #
        self.rbc2gas, self.rbc2gas_binning = disBinning(discomp       = self.rbc,
                                                        gas_highSNR   = abs(self.gas_highSNR),
                                                        bin_threshold = thre_rbc_rat,
                                                        mask          = self.mask_reg_vent,
                                                        cor           = cor_TE90*cor_flipoff)
        # from GX_defineColormaps import thre_vent
        # from GX_Map_utils import gasBinning
        # rbc2gas = np.zeros(np.shape(self.rbc))
        # mask = self.mask_reg_vent
        # rbc2gas[mask] = np.divide(self.rbc[mask],abs(self.gas_highSNR[mask]))
        # rbc2gas[rbc2gas<0] = 1e-5
        # self.rbc2gas = rbc2gas*cor_TE90*cor_flipoff
        #
        # _, self.rbc2gas_binning,_ = gasBinning(gas_highreso  = rbc2gas,
        #                                         bin_threshold = thre_vent,
        #                                         mask          = self.mask_reg_vent,
        #                                         percentile    = 99)

    def generateStats(self):
        ## calculate statistics
        from GX_Map_utils import binStats

        gas_stats = binStats(rawdata = self.ventilation,
                             bindata = self.ven_binning,
                             mask = self.mask_reg,
                             mask_all=self.mask_reg,
                             key = 'ven',
                             recondata=abs(self.gas_highreso)) # for SNR

        # the percentage is calculated by all mask, while the mean is calculated with the vent mask
        barrier_stats = binStats(rawdata = self.bar2gas,
                                 bindata = self.bar2gas_binning,
                                 mask = self.mask_reg_vent,
                                 mask_all = self.mask_reg,
                                 key = 'bar',
                                 recondata=self.barrier)

        # the percentage is calculated by all mask, while the mean is calculated with the vent mask
        rbc_stats = binStats(rawdata = self.rbc2gas,
                             bindata = self.rbc2gas_binning,
                             mask = self.mask_reg_vent,
                             mask_all = self.mask_reg,
                             key = 'rbc',
                             recondata=self.rbc)

        self.stats_box = gas_stats
        self.stats_box.update(barrier_stats)
        self.stats_box.update(rbc_stats)
        # the unit of inflation is cc
        self.stats_box['inflation'] = np.multiply(np.sum(self.mask_reg),self.FOV**3/np.shape(self.mask_reg)[0]**3)

    def generateFigures(self):
        ## make montage, plot histogram, and generate report

        from GX_Map_utils import makeMontage, makeHistogram, decideStartInterval

        from GX_defineColormaps import short_index2color, long_index2color, venhistogram_rat, barhistogram_rat, rbchistogram_rat

        ind_start, ind_inter = decideStartInterval(mask = self.mask_reg)

        ## make montage
        ven_montage = makeMontage(bin_index = self.ven_binning,
                                  ute_reg  = self.ute_reg,
                                  index2color = short_index2color,
                                  ind_start = ind_start,
                                  ind_inter = ind_inter,
                                  Subject_ID = self.Subject_ID,
                                  mon_name = self.data_dir+'/ven_montage.png')

        bar_montage = makeMontage(bin_index = self.bar2gas_binning,
                                  ute_reg = self.ute_reg,
                                  index2color = long_index2color,
                                  ind_start = ind_start,
                                  ind_inter = ind_inter,
                                  Subject_ID = self.Subject_ID,
                                  mon_name = self.data_dir+'/bar_montage.png')

        rbc_montage = makeMontage(bin_index = self.rbc2gas_binning,
                                  ute_reg = self.ute_reg,
                                  index2color = short_index2color,
                                  ind_start = ind_start,
                                  ind_inter = ind_inter,
                                  Subject_ID = self.Subject_ID,
                                  mon_name = self.data_dir+'/rbc_montage.png')

        ## make histogram
        venhistogram_rat['data'] = self.ventilation[self.mask_reg]
        venhistogram_rat['hist_name'] = self.data_dir+'/ven_hist.png'
        makeHistogram(**venhistogram_rat)

        barhistogram_rat['data'] = self.bar2gas[self.mask_reg_vent]
        barhistogram_rat['hist_name'] = self.data_dir+'/bar_hist.png'
        makeHistogram(**barhistogram_rat)

        rbchistogram_rat['data'] = self.rbc2gas[self.mask_reg_vent]
        rbchistogram_rat['hist_name'] = self.data_dir+'/rbc_hist.png'
        makeHistogram(**rbchistogram_rat)

        # barhistogram['data'] = self.bar2gas[self.mask_reg_vent]
        # barhistogram['hist_name'] = self.data_dir+'/bar_hist.png'
        # makeHistogram(**barhistogram)
        #
        # rbchistogram['data'] = self.rbc2gas[self.mask_reg_vent]
        # rbchistogram['hist_name'] = self.data_dir+'/rbc_hist.png'
        # makeHistogram(**rbchistogram)

    def generateHtmlPdf(self):

        from GX_Map_utils import genHtmlPdf

        genHtmlPdf(Subject_ID = self.Subject_ID,
                   data_dir = self.data_dir,
                   RBC2barrier = self.RBC2barrier,
                   stats_box = self.stats_box)

    def saveMat(self):

        sio.savemat(self.data_dir+'/'+self.Subject_ID+'_GX.mat',vars(self))


if __name__ == "__main__":

    if (len(sys.argv) == 2):
        data_dir = '/media/sharedrive/shared/team_documents/2018_Preclinical_PAH//RawData/HEALTHY/'+sys.argv[1]
        Subject_ID = sys.argv[1]

    elif(len(sys.argv) == 3):
        data_dir = sys.argv[1]
        Subject_ID = sys.argv[2]

    else:
        print "Usage 1: python GX_classmap.py <data directory/Subject_ID>"
        print "Usage 2: python GX_classmap.py <data directory> <Subject_ID>"
        sys.exit(-1)

    # Create helper object
    from GX_Map_utils import fullMontage
    rat = GXRat(data_dir=data_dir, Subject_ID=Subject_ID)
    rat.saveMat()
    pdb.set_trace()

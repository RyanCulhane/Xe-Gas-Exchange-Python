import numpy as np
import numpy as np
import nibabel as nib
import scipy.io as sio
import glob, os
import sys
import time

from matplotlib import pyplot as plt

import pdb

class GXSubject(object):
    'A gas exchange subject'

    def __init__(self, data_dir, Subject_ID):
        print("*********************Initiate subject. Today is a beautiful day!")
        self.data_dir = data_dir
        self.RBC2barrier = 0.0
        self.Subject_ID = Subject_ID
        self.FOV = 40.0
        self.TE90 = 460 # will be replaced later from Dixon

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

        self.stats_box = {}

    def GXRecon(self):

        print("*********************Read in 129Xe data")

        self.spectFit()
        self.reconXe()
        # self.readinXe()

        print("*********************Mask Segmentation")
        self.uteSegmentation()
        self.alignImages()
        # self.checkAlignment()
        self.uteRegister()

    def GXMapping(self):
        print("*********************Gas BiasFied Correction, Binning and Ventilation Mask")
        self.gasBiasFieldCor()
        self.gasBinning()

        print("*********************Dixon and binning")
        self.dixonDecomp()
        self.barBinning()
        self.rbcBinning()

        self.generateStats()

        print("*********************Clnical Report")
        self.generateFigures()
        self.generateHtmlPdfPPT()

        print('*********************finished program')

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

    def reconXe(self):

        from GX_Recon_utils import recon_ute, recon_dixon

        ute_file = glob.glob(self.data_dir+'/meas*1H_BHUTE_Radial*.dat')[0]
        self.ute = abs(recon_ute(ute_file))

        dixon_file = glob.glob(self.data_dir+'/meas*Xe_Radial_Dixon*.dat')[0]
        self.gas_highSNR, self.gas_highreso, self.dissolved, self.TE90 = recon_dixon(twix_file=dixon_file)
        # pdb.set_trace()

    def readinXe(self):
        ## temporal usage
        # read in dissolved and gas Xe
        fdata = self.data_dir+'/Sub'+self.Subject_ID+'_data.mat'

        mat_input = sio.loadmat(fdata)
        self.gas_highreso = mat_input['gasVol_highreso']
        self.gas_highSNR = mat_input['gasVol_highSNR']
        self.dissolved = mat_input['dissolvedVol']
        self.ute = abs(mat_input['uteVol'])

        fmask = glob.glob(self.data_dir+'/BHUTE_Sub'+self.Subject_ID+'_FID*mask*.nii')[0]
        self.mask = np.array(nib.load(fmask).get_data())

    def spectFit(self):
        from GX_Spec_utils import spect_fit

        cali_file = glob.glob(self.data_dir+'/meas*Xe_fid_cali*.dat')[0]
        dixon_file = glob.glob(self.data_dir+'/meas*Xe_Radial_Dixon*.dat')[0]

        self.RBC2barrier, self.spectfit_box = spect_fit(twix_dixon_file=dixon_file, twix_cali_file=cali_file, Subject_ID = self.Subject_ID)

    def uteSegmentation(self):
        ## temporal usage
        from GX_CNNpredict import CNNpredict

        # fmask = 'BHUTE_Sub002102_FID49886_mask_grow.nii'
        # self.mask = np.array(nib.load(fmask).get_data(),dtype='bool')
        # fute = glob.glob(self.data_dir+'/BHUTE_Sub'+self.Subject_ID+'_FID*recon.nii')[0]
        # fute = 'BHUTE_Sub002102_FID49886_recon.nii'
        # self.ute = np.array(nib.load(fute).get_data())

        self.mask = CNNpredict(ute = self.ute)

    def alignImages(self):

        def alignrot(x):
            return(np.flip(np.flip(np.flip(x,0),1),2))

        self.ute = alignrot(self.ute)
        self.mask = alignrot(self.mask)
        self.gas_highreso = alignrot(self.gas_highreso)
        self.gas_highSNR = alignrot(self.gas_highSNR)
        self.dissolved = alignrot(self.dissolved)

    def uteRegister(self):
        ## temporal usage
        from GX_Map_utils import register

        self.ute_reg, self.mask_reg = register(gas_highreso = abs(self.gas_highreso),
                                               ute          = self.ute,
                                               mask         = self.mask.astype(float))

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
        from GX_defineColormaps import thre_vent
        from GX_Map_utils import gasBinning

        self.ventilation, self.ven_binning, self.mask_reg_vent = gasBinning(gas_highreso  = abs(self.gas_highreso_cor),
                                                                            bin_threshold = thre_vent,
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
        from GX_defineColormaps import thre_bar
        from GX_Map_utils import disBinning

        cor_TE90 = np.exp(self.TE90/2000.0)/np.exp(self.TE90/50000.0)
        cor_flipoff = 100*np.sin(0.5*np.pi/180)/np.sin(20*np.pi/180)

        self.bar2gas, self.bar2gas_binning = disBinning(discomp       = self.barrier,
                                                        gas_highSNR   = abs(self.gas_highSNR),
                                                        bin_threshold = thre_bar,
                                                        mask          = self.mask_reg_vent,
                                                        cor           = cor_TE90*cor_flipoff)

    def rbcBinning(self):
        ## binning for barrier
        from GX_defineColormaps import thre_rbc
        from GX_Map_utils import disBinning

        cor_TE90 = np.exp(self.TE90/2000.0)/np.exp(self.TE90/50000.0)
        cor_flipoff = 100*np.sin(0.5*np.pi/180)/np.sin(20*np.pi/180)

        self.rbc2gas, self.rbc2gas_binning = disBinning(discomp       = self.rbc,
                                                        gas_highSNR   = abs(self.gas_highSNR),
                                                        bin_threshold = thre_rbc,
                                                        mask          = self.mask_reg_vent,
                                                        cor           = cor_TE90*cor_flipoff)

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
        self.stats_box['inflation'] = np.multiply(np.sum(self.mask_reg),self.FOV**3/np.shape(self.mask_reg)[0]**3/1000.0)

    def generateFigures(self):
        ## make montage, plot histogram, and generate report

        from GX_Map_utils import makeMontage, makeHistogram, decideStartInterval

        from GX_defineColormaps import short_index2color, long_index2color, venhistogram, barhistogram, rbchistogram

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
        venhistogram['data'] = self.ventilation[self.mask_reg]
        venhistogram['hist_name'] = self.data_dir+'/ven_hist.png'
        makeHistogram(**venhistogram)

        barhistogram['data'] = self.bar2gas[self.mask_reg_vent]
        barhistogram['hist_name'] = self.data_dir+'/bar_hist.png'
        makeHistogram(**barhistogram)

        rbchistogram['data'] = self.rbc2gas[self.mask_reg_vent]
        rbchistogram['hist_name'] = self.data_dir+'/rbc_hist.png'
        makeHistogram(**rbchistogram)

    def generateHtmlPdfPPT(self):

        from GX_Map_utils import genHtmlPdf
        from GX_defineColormaps import referece_stats

        genHtmlPdf(Subject_ID = self.Subject_ID,
                   data_dir = self.data_dir,
                   RBC2barrier = self.RBC2barrier,
                   stats_box = self.stats_box,
                   referece_stats = referece_stats)

    def saveMat(self):

        sio.savemat(self.data_dir+'/'+self.Subject_ID+'.mat',vars(self))

    def sendEmail(self):
        from GX_email import send_email

        send_email(data_dir=self.data_dir, Subject_ID=self.Subject_ID)

    def readFromMat(self):
        # read in all members from a saved matlab file
        fdata = self.data_dir+'/'+self.Subject_ID+'.mat'

        if not os.path.isfile(fdata):
            raise Exception('Mat file of this subject does not exist, need to create it')

        mat_input = sio.loadmat(fdata)

        self.RBC2barrier = mat_input['RBC2barrier'].flatten()
        self.TE90 = mat_input['TE90'].flatten() # will be replaced later from Dixon

        self.gas_highreso = mat_input['gas_highreso']
        self.gas_highSNR = mat_input['gas_highSNR']
        self.gas_highreso_cor = mat_input['gas_highreso_cor']
        self.gas_biasfield = mat_input['gas_biasfield']
        self.dissolved = mat_input['dissolved']

        self.rbc = mat_input['rbc']
        self.barrier = mat_input['barrier']

        self.ventilation = mat_input['ventilation']
        self.rbc2gas = mat_input['rbc2gas']
        self.bar2gas = mat_input['bar2gas']
        self.ven_binning = mat_input['ven_binning']
        self.rbc2gas_binning = mat_input['rbc2gas_binning']
        self.bar2gas_binning = mat_input['bar2gas_binning']

        self.ute = mat_input['ute']
        self.mask = mat_input['mask'].astype(bool)
        self.ute_reg = mat_input['ute_reg']
        self.mask_reg = mat_input['mask_reg'].astype(bool)
        self.mask_reg_vent = mat_input['mask_reg_vent'].astype(bool)

        self.stats_box = mat_input['stats_box']

if __name__ == "__main__":

    if (len(sys.argv) == 2):
        data_dir = '/home/ziyiw/Patients/'+sys.argv[1]
        Subject_ID = sys.argv[1]

    elif(len(sys.argv) == 3):
        data_dir = sys.argv[1]
        Subject_ID = sys.argv[2]

    else:
        print "Usage 1: python GX_classmap.py <data directory/Subject_ID>"
        print "Usage 2: python GX_classmap.py <data directory> <Subject_ID>"
        sys.exit(-1)

    # pdb.set_trace()
    # Create helper object
    # from GX_Map_utils import fullMontage
    start = time.time()
    subject = GXSubject(data_dir=data_dir, Subject_ID=Subject_ID)
    # subject.readFromMat()
    subject.GXRecon()
    subject.GXMapping()
    subject.saveMat()
    endt = time.time()
    print(endt-start)

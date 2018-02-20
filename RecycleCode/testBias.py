import SimpleITK as sitk
import numpy as np
import nibabel as nib
import scipy.io as sio
import pdb
from GX_utils import fullMontage

fdata = 'Sub002102_data.mat'

mat_input = sio.loadmat(fdata)
gas_highreso = abs(mat_input['gasVol_highreso'])

fmask = 'BHUTE_Sub002102_FID49886_mask_grow.nii'
mask = np.array(nib.load(fmask).get_data())

sitk_mask = sitk.GetImageFromArray(mask)
sitk_image = sitk.GetImageFromArray(gas_highreso)

bth = sitk.BinaryThresholdImageFilter()
sitk_mask_b = bth.Execute(sitk_mask)

corrector = sitk.N4BiasFieldCorrectionImageFilter()

# sitk_image = sitk.Cast( sitk_image, sitk.sitkFloat32 )
# corrector.SetMaximumNumberOfIterations()

pdb.set_trace()
sikt_image_cor = corrector.Execute( sitk_image, sitk_mask_b )

# image_cor = sitk.GetArrayFromImage(sitk_image_cor)

pdb.set_trace()

import os
import pdb
import numpy as np
import nibabel as nib

def biasFieldCor(image, mask):

    pathInput = 'image.nii'
    pathMask = 'mask.nii'
    pathOutput = 'image_cor.nii'
    pathBiasField = 'biasfield.nii'

    pathN4 = './N4BiasFieldCorrection'

    nib.save(abs(image), pathInput)
    nib.save(mask, pathMask)
    
    # cmd = pathN4+' -d 3 -i '+pathInput+' -s 2 -x '+pathMask+' -o ['+pathOutput+', '+pathBiasField+']'
    cmd = pathN4+' -d 3 -i '+pathInput+' -x '+pathMask+' -o ['+pathOutput+', '+pathBiasField+']'

    os.system(cmd)

    image_cor = np.array(nib.load(pathOutput).get_data())
    image_biasfield = np.array(nib.load(pathBiasField).get_data())

    os.remove(pathOutput)
    os.remove(pathBiasField)

return image_cor, image_biasfield

"""
nifti_io.py
Used for NIfTI-1 file IO
"""
import os
import sys
import nibabel as nib
import numpy as np
import pdb

# TODO: look into the documentation and find a way to preserve meta data
def write_nii(array_data, filename, path = "", affine = None):
    """write np array into nii file"""
    if affine is None:
        print("No information about the global coordinate system")
        affine = np.diag([1,1,1,1])
        filename = "No_GCS_" + filename
    #pdb.set_trace()
    array_data = np.int16(array_data)
    array_img = nib.Nifti1Image(array_data, affine)
    save_fid = os.path.join(path,filename)
    array_img.to_filename(filename)
#    pdb.set_trace()
    return filename

def read_nii_image(input_fid):
    """read the nii image data into numpy array"""
    img = nib.load(input_fid)
    return img.get_data()

def read_nii_object(input_fid):
    """ directly read the nii object """
    return nib.load(input_fid)

# write resampled nii file. Keep the corresponding metadata unchanged, for 3dslicer visualization
def write_resampled_nii(original_obj, data_vol, filename, new_res = None,debug = False):
    """ write resampleed nii file """
    affine = original_obj.get_affine()
    
    if debug == True:
        old_affine = _affine
    # replace the affine scale parameter with new resolution
    if new_res != None:
        for i in range(3):
            #pdb.set_trace()
            affine[i,i] = new_res[i]
        affine[0,0] *= -1
        
    if debug == True:
        print("Old affine matrix: ")
        print(old_affine)
        print("New affine matrix: ")
        print(affine)
    
    output_obj = nib.Nifti1Image(data_vol, affine)
    output_obj.to_filename("resampled_" + filename)
    return "resampled_" + filename
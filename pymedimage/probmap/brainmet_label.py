"""
Generate a binary mask given the RTStruct file and perform gaussian filtering for brain_met
probability map generation

This should be done in after dicom files are 
"""
# TODO: organize previously written nifti operation code (brainMets) and this one
import pymedimage.dcmio as dio
import numpy as np
import scipy.ndimage.filters as sfilter
import nibabel as nib
import showvolume as viz
import os
import dicom as dcm 
import pymedimage.rttypes as rts
import copy
# TODO: unify the usage of BaseVolume and dense np array
# Use the same base data structure BaseVolume to store masked 
# TODO: Consider one thing: the dataset is stored in nii files, while the 
# rtstruct is in dicom
class RTContourLabel(rts.ROI):
""" Generating training labels from rtstruct. Supporting pixel-wise dense label map generation and soft probability
    map generation given the location of the contour center.

    extenstion of pymedimage.rttypes.ROI class  
    Member objects:
        self.reference_slice_fid: fid of the dicom image slice which it attaches to
        self.contour_loc_idx: index of the center of contour center in self.array
        self.contour_loc_coord: coordinate of the center of contour center
        self.reference_slice_FOR: rttypes.FrameOfReference of the target volume where the rtstruct is resampled to 
        self.array 
    Methods:
        self._find_center()
        self._loc_to_dense_array() : generate a location point array with the same size of reference slice
        self._coord2idx() : convert dicom coordinate to numpy index given FrameOfReference
        self._idx2coord() : convert numpy index to dicom coordinate given FrameOfReference 
        self.smooth_prob_map() : return the smoothed probability map with the same size of reference dicom
        self.tumer_mask(): return pixelwise tumer mask with the same size of reference dicom 
        
"""

    def __init__(self, base_ROI, reference_slice_fid = None):
        """ initialize from a pymedimage.rttypes.BaseVolume object """
        self.__dict__ = copy.deepcopy(base_ROI.__dict__)
        self.contour_loc_idx = self._find_ctr()
        self.contour_loc_coord = self._idx2coord(self.contour_loc_idx, self.frameofreference)
        self.array = super(RTContourLabel, self).makeDenseMask()
        if reference_slice_fid is not None:
            self.reference_slice_fid = reference_slice_fid 
       
    def _find_center(self):
        """ find the center of this contour 
            Assuming that backgound is 0 while foreground is not
        """
        print("Warning: The ROI center calculation code cannot deal with the case where multiple connected component are present")
        _foreground_idx = np.where(self.array > 0)
        return [ (np.mean(_foreground_idx[0])), (np.mean(_foreground_idx[1])), (np.mean(_foreground_idx[2])) ]

    def _center_array(self):
        idx = tuple([int(axis) for axis in self.contour_loc])
        _ctr_array = np.zeros(self.array.shape)
        _ctr_array[idx] = 1
        return _ctr_array

    def _coord2idx(self, coord, reference):
        """ conversion between dicom coordinate and numpy indices
            Args:
                reference: FrameOfReference object
        """
        np_dims = []
        for idx in range(len(reference.spacing)):
            np_dims.append( int(( coord[idx] - reference.starting[idx] ) / reference.spacing[idx]) )
        return tuple(reversed(np_dims))

    def _idx2coord(self, idx, reference):
        """ reverse of _coord2idx"""
        dcm_coords = []
        for axis in range(len(reference.spacing)):
            dcm_coords.append( float(idx[axis] * reference.spacing[axis] + ref.starting[axis]) )
        return tuple(reversed(dcm_coords))

    def _loc_to_dense_array(self, ref_fid = None):
        """ given the reference slice, generate a mask with the same size """
        if ref_fid is not None:
            _ref_filename = ref_fid
        else:
            _ref_filename = self.reference_slice_fid
        # assume that all the other slices of this volume are saved in a same directory
        if self.reference_slice_FOR is None:
            self.reference_slice_FOR = rts.FrameOfReference.from_dcm_fid(_ref_filename)
        # Note: This is quite beautiful since similar implementation has been also used in the parent class
        # Note: The order of numpy array axis and dicom are opposite
        _npdim2, _npdim1, _npdim0 = self.reference_slice_FOR.size
        # convert the center coordinate into index in the dense volume 
        dense_location_map = np.zeros([_npdim0, _npdim1, _npdim2])
        dense_contour_center_idx = self._coord2idx(self.contour_loc_coord, self.reference_slice_FOR) 
        dense_location[dense_contour_center_idx] = 1
        return dense_location_map

    # TODO: Add support of other kernels
    def smooth_prob_map(self, sigma = None, ref_fid = None, kernel = 'Gaussian'):
        """ Generate a soft probability map of contour occurance 
            Args:
                Sigma: the gaussian smoothing radias in pixels
                ref_fid: dicom file fid of image volume where the contour will be conformed to 
                kernel: certain kernel to be used other than gaussian
        """
        dense_location_map = self._loc_to_dense_array(ref_fid)
        if kernel = 'Gaussian':
            if sigma = None:
                sigma = 5  
            return sfilter.guassian_filter(dense_location_map, sigma, mode = 'nearest')
        else:
            raise Exception(" This function is still under construction! ")        

        


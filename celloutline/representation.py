# encoding: utf-8
""" Representation stores and compares cell representations 
Author: CDW
"""

# Standard or installed
import numpy as np
import scipy.ndimage
# Local
from . import conversions  # does heavy lifting of binary->other->back


__all__ = ["BinaryVoxel", "SpreadVoxel",
           "SpiralizedTrace", "mesh_error"]


""" Comparisons between the representations; happen at mesh level """


def mesh_error(mesh1, mesh2):
    """Error (intersection over union) of the two meshes"""
    intersection = mesh1.intersection(mesh2)
    union = mesh1.union(mesh2)
    mesh_error = intersection.volume/union.volume
    return mesh_error


""" Representational classes """
 

class Representation:
    """What all the forms of a segmentation are based on"""
    def __init__(self, id, source=None):
        """Remember where we came from
        
        Parameters
        ----------
        id: string
            key or filename
        source: collection or None
            dict or other collection for key, None if the id 
            is a filename
        """
        self.id = id
        self.source = source
        self._voxels = None
        self._spread = None
        self._spiral = None
        self._mesh = None
    
    def error(self, other):
        """Mismatch between this and another representation"""
        return mesh_error(self.mesh, other.mesh)


class BinaryVoxel(Representation):
    """Reading and conversion of binary cell segmentations"""
    def __init__(self, id, source=None, voxels=None):
        """This is the native representation of 3D microscopy data

        We read in the binary voxel representation from numpy arrays 
        that have been created from the original microscopy data by 
        alignment and downsampling to uniform square voxels.

        Parameters
        ----------
        id: string
            key or filename
        source: collection or None
            dict or other collection for key, None if the id 
            is a filename
        voxels: i-by-j-by-k array or None
            the voxels of interest or nothing (in which case we load
            them on request)
        """
        super().__init__(id, source)
        self._voxels = voxels # lazy load initial data

    @property
    def voxels(self):
        """The binary voxel representation of the cell"""
        if self._voxels is None:
            if self.source is None: # load from disc
                self._voxels = np.load(self.id)['cell']
            else: # load from source collection
                self._voxels = self.source[self.id]
        return self._voxels

    @property
    def mesh(self):
        """The mesh of the voxels"""
        if self._mesh is None:
            self._mesh = conversions.binary_to_mesh(self.voxels)
            assert self._mesh.is_watertight, "leaky mesh"
        return self._mesh

    @property
    def spread(self):
        """A spread representation of the binary voxels"""
        if self._spread is None:
            spread_voxels = conversions.binary_to_spread(self.voxels)
            self._spread = SpreadVoxels(self.id, 
                                        self.source, 
                                        spread_voxels)
        return self._spread

    def spiral(self, unitspiral=None, num_pts=500):
        """A spiral representation of the cell"""
        if self._spiral is None:
            spidict = conversions.binary_to_spiral(
                self.voxels, unitspiral, num_pts)
            self._spiral = SpiralizedTrace(
                self.id, self.source, **spidict)
        return self._spiral
        

class SpiralizedTrace(Representation):
    """Spiral trace along the cell surface"""
    def __init__(self, id, source, radii, unitspiral, origin):
        """ Contain a representation of a spiralized trace and 
        the unit spiral that is needed to generate it. 
        
        Parameters
        ----------
        id: string
            key or filename
        source: collection or None
            dict or other collection for key, None if the id 
            is a filename
        radii: 1-by-n array
            list of distances from the origin to the first shell 
            intersection for a given spiral
        unitspiral: UnitSpiral object
            contains unit rays (in same order as radii) with angles
        origin: 3-by-1 array
            xyz offset of the center of the segmentation
        """
        super().__init__(id, source)
        self._point_cloud = None
        self.spiral_dict = {'unitspiral': unitspiral,
                            'radii': radii,
                            'origin': origin}

    @property
    def unitspiral(self):
        return self.spiral_dict['unitspiral']

    @property
    def radii(self):
        return self.spiral_dict['radii']

    @property
    def origin(self):
        return self.spiral_dict['origin']

    @property
    def point_cloud(self):
        if self._point_cloud is None:
            self._point_cloud = conversions.spiral_to_point_cloud(
                **self.spiral_dict)
        return self._point_cloud

    @property
    def mesh(self):
        if self._mesh is None:
            self._mesh = conversions.spiral_to_mesh(**self.spiral_dict)
        return self._mesh


class SpreadVoxels(Representation):
    """Level set derived from binary voxels"""
    def __init__(self, id, source, spread, cutoff=0.0):
        """
        Parameters
        ----------
        id: string
            key or filename
        source: collection or None
            dict or other collection for key, None if the id 
            is a filename
        spread: i-by-j-by-k array
            the spread voxels of interest 
        cutoff: float
            where we set the isosurface for the mesh conversion
        """
        super().__init__(id, source)
        self._spread = spread 
        self._cutoff = cutoff
        self._cutoff_change = False

    @property
    def cutoff(self):
        return self._cutoff

    @cutoff.setter
    def cutoff(self, newcutoff):
        """Set mesh conversion cutoff and trigger new mesh"""
        self._cutoff = newcutoff
        self._cutoff_change = True

    @property
    def mesh(self):
        """New mesh if not done or if cuttoff has changed"""
        if self._mesh is None or self._cutoff_change:
            self._mesh = conversions.spread_to_mesh(
                self._spread, self._cutoff)
            self._cutoff_change = False
        return self._mesh

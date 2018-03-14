# encoding: utf-8
""" Representation stores and compares cell representations
Author: CDW
"""

# Standard or installed
import numpy as np
# Local
from . import conversion  # does heavy lifting of binary->other->back


__all__ = ("BinaryVoxel", "SpreadVoxel",
           "SpiralizedTrace", "mesh_error")


""" Comparisons between the representations; happen at mesh level """


def mesh_error(mesh1, mesh2):
    """Error (intersection over union) of the two meshes"""
    intersection = mesh1.intersection(mesh2)
    union = mesh1.union(mesh2)
    error = intersection.volume/union.volume
    return error


""" Representational classes """


class Representation:
    """What all the forms of a segmentation are based on"""
    def __init__(self, name, source=None):
        """Remember where we came from

        Parameters
        ----------
        name: string
            key or filename
        source: collection or None
            dict or other collection for key, None if the name
            is a filename
        """
        self.name = name
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
    def __init__(self, name, source=None, voxels=None):
        """This is the native representation of 3D microscopy data

        We read in the binary voxel representation from numpy arrays
        that have been created from the original microscopy data by
        alignment and downsampling to uniform square voxels.

        Parameters
        ----------
        name: string
            key or filename
        source: collection or None
            dict or other collection for key, None if the name
            is a filename
        voxels: i-by-j-by-k array or None
            the voxels of interest or nothing (in which case we load
            them on request)
        """
        super().__init__(name, source)
        self._voxels = voxels  # lazy load initial data

    @property
    def voxels(self):
        """The binary voxel representation of the cell"""
        if self._voxels is None:
            if self.source is None:  # load from disc
                self._voxels = np.load(self.name)['cell']
            else:  # load from source collection
                self._voxels = self.source[self.name]
        return self._voxels

    def mesh(self, step=1):
        """The mesh of the voxels"""
        if self._mesh is None:
            self._mesh = conversion.binary_to_trimesh(self.voxels, step)
            assert self._mesh.is_watertight, "leaky mesh"
        return self._mesh

    def spread(self):
        """A spread representation of the binary voxels"""
        if self._spread is None:
            spread_voxels = conversion.binary_to_spread(self.voxels)
            self._spread = SpreadVoxel(self.name,
                                       self.source,
                                       spread_voxels)
        return self._spread

    def spiral(self, unitspiral=None, num_pts=500):
        """A spiral representation of the cell"""
        if self._spiral is None or unitspiral is not None:
            spidict = conversion.binary_to_spiral(
                self.voxels, unitspiral, num_pts)
            self._spiral = SpiralizedTrace(
                self.name, self.source, **spidict)
        return self._spiral

    def to_vector(self):
        """Row vector representation for use with models"""
        return self.voxels.flatten()

    def from_vector(self, vector):
        """Reshape in an input vector back to a voxel array.
        Return a new BinaryVoxel instance for the new data
        """
        newvox = vector.reshape(self.voxels.shape)
        return BinaryVoxel(self.name, 'reconstruction', newvox)


class SpiralizedTrace(Representation):
    """Spiral trace along the cell surface"""
    def __init__(self, name, source, radii, unitspiral, origin):
        """ Contain a representation of a spiralized trace and
        the unit spiral that is needed to generate it.

        Parameters
        ----------
        name: string
            key or filename
        source: collection or None
            dict or other collection for key, None if the name
            is a filename
        radii: 1-by-n array
            list of distances from the origin to the first shell
            intersection for a given spiral
        unitspiral: UnitSpiral object
            contains unit rays (in same order as radii) with angles
        origin: 3-by-1 array
            xyz offset of the center of the segmentation
        """
        super().__init__(name, source)
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
            self._point_cloud = conversion.spiral_to_point_cloud(
                **self.spiral_dict)
        return self._point_cloud

    @property
    def mesh(self):
        if self._mesh is None:
            self._mesh = conversion.spiral_to_trimesh(**self.spiral_dict)
        return self._mesh

    def to_vector(self):
        """Row vector representation for use with models"""
        return self.radii

    def from_vector(self, vector):
        """Reshape in an input vector back to a voxel array.
        Return a new SpiralizedTrace instance for the new data
        """
        new_spiral_dict = self.spiral_dict.copy()
        new_spiral_dict['radii'] = vector
        return SpiralizedTrace(self.name, 'reconstruction', **new_spiral_dict)


class SpreadVoxel(Representation):
    """Level set derived from binary voxels"""
    def __init__(self, name, source, spread, cutoff=0.0):
        """
        Parameters
        ----------
        name: string
            key or filename
        source: collection or None
            dict or other collection for key, None if the name
            is a filename
        spread: i-by-j-by-k array
            the spread voxels of interest
        cutoff: float
            where we set the isosurface for the mesh conversion
        """
        super().__init__(name, source)
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
            self._mesh = conversion.spread_to_trimesh(
                self._spread, self._cutoff)
            self._cutoff_change = False
        return self._mesh

    def to_vector(self):
        """Row vector representation for use with models"""
        return self._spread.flatten()

    def from_vector(self, vector):
        """Reshape in an input vector back to a voxel array.
        Return a new BinaryVoxel instance for the new data
        """
        newspread = vector.reshape(self._spread.shape)
        return BinaryVoxel(self.name, 'reconstruction', newspread)

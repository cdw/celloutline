#!/usr/bin/env python
# encoding: utf-8
""" Representation converts between and compares cell representations 
Author: CDW
"""

# Standard or installed
import numpy as np
import trimesh #handles mesh creation and interface to OpenSCAD
from skimage.measure import marching_cubes
import scipy.ndimage
import scipy.spatial
# Local
from .geom import sphere_to_cart, sprial_to_point_cloud, Spiral

__all__ = ["BinaryVoxel", "SpreadVoxel", "SpiralizedTrace"]

""" Raw representation conversions, one type of encoding to the next

We want to support the following transitions:
 - binary -> (mesh, spiral, spread)
 - spiral -> (point-cloud, mesh)
 - spread -> (binary, mesh)
"""

def binary_to_mesh(binary):
    """ Convert binary voxel image to a mesh

    Marching cubes meshed output from 
    http://scikit-image.org/docs/dev/api/skimage.measure.html#marching-cubes-lewiner
    """
    verts, faces, _, _ = marching_cubes(binary)
    tm = trimesh.Trimesh(verts, faces)
    return Mesh(tm)

def binary_to_spiral(binary, spiral=None, num_pts=500):
    """Convert a binary voxel image into a cell spiral
    
    Parameters
    ---------
    binary: i-by-j-by-k array
        binary voxel segmentation
    spiral: None or Spiral object
        spiral class to provide rays for intersections
    num_pts: None or int (default 500)
        number of points to use if we need to create our own spiral

    Returns
    -------
    xyz: num_pts-by-3 array
        list of xyz values for cell spiral intersections
    radii: num_pts-by-3
        distance to each output point from the origin
    origin: 1-by-3 array
        xyz origin of the spiral
    """
    if spiral is None:
        spiral = Spiral(num_pts)
    # Find shell to intersect with
    shell = segmentation.find_boundaries(binary, connectivity=1, mode='outer')
    surf_coords = list(np.array(shell.nonzero()).T)
    origin = np.divide(shell.shape,2).astype(int)
    # Find intersections
    xyz = []
    for ray in spiral.xyz:
        xyz.append(nearest_intersecting((origin, ray), surf_coords))
    radii = [dist(origin, pt, 0.5) for pt in output_pts]
    return xyz, radii, origin

def binary_to_spread(binary):
    """Convert a binary voxel image into a levelset
    positive inside cell and neg outside
    """
    # The interior is dist to shell
    pos_interior = scipy.ndimage.distance_transform_edt(binary)
    # The exterior is -1 * dist to shell
    neg_exterior = -scipy.ndimage.distance_transform_edt(1-binary)
    # Combine to the levelset
    spread = pos_interior + neg_exterior
    return spread

def spiral_to_point_cloud(radii, Spiral, origin):
    """Convert a spiral trace back to a xyz point cloud given a set of rays"""
    # Like the original spiral, but with new radii
    rpt = np.hstack((np.expand_dims(radiim -1), Spiral.rpt[:,1:]))
    # Convert the whole thing to cartesian and offset by the origin
    xyz = sphere_to_cart(rpt)
    xyz += origin
    return xyz

def spiral_to_mesh(radii, Spiral, origin):
    """Convert a spiral to a mesh via the mesh of the spiral's rays"""
    # Convert radii to xyz points
    point_cloud = spiral_to_point_cloud(radii, Spiral, origin)
    # How would we mesh the original unit-circle spiral?
    faces = scipy.spatial.ConvexHull(Spiral.rpt).simplices
    # OK, well do that again
    mesh = trimesh.Trimesh(point_cloud, faces)
    # Flip the faces to be all pointing outwards
    mesh.fix_normals()
    return mesh

def spread_to_mesh(spread, cutoff):
    """Convert levelset to mesh by way of binary"""
    return binary_to_mesh(spread_to_binary(spread, cutoff))


""" Comparisons between the representations

All comparisons happen at mesh level, the paths to mesh are
 - binary -> mesh
 - binary -> spiral -> point cloud -> mesh
 - binary -> spread -> mesh
"""

def mesh_mismatch(mesh1, mesh2):
        """The intersection over the union of the two meshes"""
    mesh1.






class Representation:
    """Do we want a base class?"""
    def __init__(self):
        return



class BinaryVoxel:
    """Reading and comparison of binary cell segmentations"""
    def __init__(self, input_fn, voxels=None):
        """This is the native representation of 3D microscopy data

        We read in the binary voxel representation from numpy arrays 
        that have been created from the original microscopy data by 
        alignment and downsampling to uniform square voxels.

        Parameters
        ----------
        input_fn: string
            path to underlying voxel file
        """
        self.fn = input_fn 
        self._voxels = voxels # lazy load initial data
        self._spread = None   # and don't create downstream
        self._spiral = None   # data until it is needed for 
        self._mesh = None     # reduced memory consumption
        # Cache results when possible, keep track
        vox_recalc = voxels is None # needed if not passed
        self.needs_recalc = {'voxels': vox_recalc,
                             'spread': True, 
                             'spiral': True,
                             'mesh': True}
        return 

    @property
    def voxels(self):
        """The binary voxel representation of the cell"""
        if self.needs_recalc['voxels'] is True:
            self._voxels = np.load(self.fn)['cell']
            self.needs_recalc['voxels'] = False
        return self._voxels

    @property
    def spread(self, options={}):
        """A spread representation of the binary voxels"""
        if 
        sl
        return SpreadVoxel

    def mismatch(self, recon):
        orig = self.read()
        size = orig.size
        assert orig.shape==recon.shape, "shapes don't match"
        absdiff = np.sum(np.abs(recon-orig))
        normed = absdiff/size
        return normed

    def spiral(self, spiral, origin, options={}):
        """ Convert binary image to radius, phi, theta spiral

        Parameters
        ----------
        spiral: Spiral class
            spiral.rpt are the rays we will trace along to find 
            intersections with the shell of the segmentation
        origin: 3 length vector
            x,y,z origin of the segmentation from which the rays emerge
        options: dict, optional
            
        Returns
        -------
        spiral: SpiralizedTrace class instance
            instance of a spiralized trace

        """
        shell = segmentation.find_boundaries(cell, connectivity=1, mode='outer')
        surf_coords = list(np.array(shell.nonzero()).T)
        origin = np.divide(shell.shape,2).astype(int)
        output_pts = []
        for ray in sample_pts:
            output_pts.append(nearest_intersecting((origin, ray), surf_coords))
        radii = [dist(origin, pt, 0.5) for pt in output_pts]
        return output_pts, radii
        
    def mesh(self, options={}):
        """ Convert binary image to a mesh

        Parameters
        ----------
        options: dict, optional
            Marching cubes control options, documented in skimage_. 
            Notable settings include, ``level``, a float level to 
            draw the boundary at.
        
        Returns
        -------
        mesh: Mesh class instance
            Marching cubes meshed output

        .. _skimage: http://scikit-image.org/docs/dev/api/skimage.measure.html#marching-cubes-lewiner
        """
        verts, faces, _, _ = marching_cubes(self.binary, **options)
        tm = trimesh.Trimesh(verts, faces)
        return Mesh(tm)



class SpiralizedTrace:
    """Reading and comparison of cell surface traces"""
    def __init__(self, inobj, rays_rpt, origin):
        """inobj must be an array or support a read method"""
        self.inobj = inobj
        self.rays_rpt = rays_rpt
        self.origin = origin

    @property
    def trace(self):
        if type(self.inobj)==np.array:
            return self.inobj
        else:
            return self.inobj.read()

    @property
    def point_cloud(self):
        # rpt = np.hstack((np.expand_dims(self.trace, -1), self.rays_rpt[:,1:]))
        # xyz = sphere_to_cart(rpt)
        # xyz += self.origin
        #return xyz
        return sprial_to_point_cloud(self.trace, self.rays_rpt, self.origin)

    def mismatch(self, shell):
        shell_locs = np.array(np.nonzero(shell)).T
        pt_cloud = self.point_cloud
        dists = [np.min(spdist.cdist(pt[np.newaxis,:], shell_locs)) for pt in pt_cloud]
        return np.mean(dists)
        pass



class SpreadVoxels:
    """SpreadVoxels is a derived structure from a set of binary voxels"""
    def __init__(self, ):
        self.parent = parent_binary_voxels
        self.needs_recalc = {'spread': True,
                             'mesh': True}
        return

    @property
    def spread(self):
        """Calculate the level-set for the parent binary voxels"""
        if self.needs_recalc['spread']:
            vox = self.parent.voxels
            neg_interior = -scipy.ndimage.distance_transform_edt(vox)
            pos_exterior = scipy.ndimage.distance_transform_edt(1-vox)
            self._spread = pos_exterior + neg_interior
            self.needs_recalc['spread'] = False
    
    def mesh_error(self, mesh, 
    @property
    def mesh:
    def read_and_transform(self, fn):
        pass
    def _write(self, output_dir):
        pass
    @staticmethod
    def mismatch(original, recon):
        pass
    def mismatch_all(self, recons):
        pass


class Mesh:
    def __init__(self, input_tm):
        self.mesh = input_tm

    def union(self, tm):
        return trimesh.union(self.mesh, tm)

    def intersection(self, tm):
        return trimesh.intersection(self.mesh, tm)

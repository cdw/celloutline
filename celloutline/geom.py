# encoding: utf-8
""" Geometric transforms and supporting concepts: consequences of 3D world
Author: CDW
"""

# Standard or installed
import numpy as np
import scipy.spatial
from numba import jit
# Local
from . import greedy

""" Coordinate conversion: xyz to rpt and back """

def cart_to_sphere(xyz): 
    """Take in xyz row vectors and return rpt row vectors
    Convention notation: rpt is radial, polar, azimuthal order 
    """
    sphere = np.zeros(xyz.shape)
    xy = xyz[:,0]**2 + xyz[:,1]**2
    sphere[:,0] = np.sqrt(xy + xyz[:,2]**2) #radial
    sphere[:,1] = np.arctan2(np.sqrt(xy), xyz[:,2]) # polar elevation defined from Z-axis down
    sphere[:,2] = np.arctan2(xyz[:,1], xyz[:,0]) # azimuthal
    return sphere


def sphere_to_cart(rpt):
    """Take in rpt row vectors and return xyz row vectors
    Convention notation: rpt is radial, polar, azimuthal order 
    """
    cart = np.zeros(rpt.shape)
    cart[:,0] = rpt[:,0] * np.cos(rpt[:,2]) * np.sin(rpt[:,1])
    cart[:,1] = rpt[:,0] * np.sin(rpt[:,2]) * np.sin(rpt[:,1])
    cart[:,2] = rpt[:,0] * np.cos(rpt[:,1])
    return cart


""" Ray intersection: support binary segmentation -> spiral """

@jit
def _intersect(ray, box):
    """Does that ray hit that box?
    
    Parameters
    ----------
    ray: (tuple of length 3, tuple of length 3)
        Ray of form ((xyz origin), (xyz unit vector))
    box: (tuple of length 3, tuple of length 3) 
        Box of form ((xyz left bottom),(xyz right top))
    
    Returns
    -------
    intersect: boolean
        True if the ray intersects the box
    """
    ray_origin, ray_direction = ray
    i_x, i_y, i_z = ray_origin
    d_x, d_y, d_z = ray_direction
    box_lb, box_rt = box
    lb_x, lb_y, lb_z = box_lb
    rt_x, rt_y, rt_z = box_rt
    # inverse of ray directions
    dirfrac_x = 1/d_x if d_x!=0 else np.inf
    dirfrac_y = 1/d_y if d_y!=0 else np.inf
    dirfrac_z = 1/d_z if d_z!=0 else np.inf
    # lb is the corner of AABB with minimal coordinates - left bottom, rt is maximal corner
    # r.org is origin of ray
    t1 = (lb_x - i_x)*dirfrac_x
    t2 = (rt_x - i_x)*dirfrac_x
    t3 = (lb_y - i_y)*dirfrac_y
    t4 = (rt_y - i_y)*dirfrac_y
    t5 = (lb_z - i_z)*dirfrac_z
    t6 = (rt_z - i_z)*dirfrac_z
    tmin = max(max(min(t1, t2), min(t3, t4)), min(t5, t6))
    tmax = min(min(max(t1, t2), max(t3, t4)), max(t5, t6))
    # if tmax < 0, ray (line) is intersecting AABB, but the whole AABB is behind us
    if tmax < 0:
        return False
    # if tmin > tmax, ray doesn't intersect AABB
    if tmin > tmax:
        return False
    return True


def dist_w_offset(origin, pt, correction=0):
    """Distance between origin and point with optional correction factor"""
    pt = np.add(pt, correction)
    distance = scipy.spatial.distance.euclidean(origin, pt)
    return distance


def nearest_intersecting(ray, voxel_corners):
    """What intersecting voxel is nearest the ray's origin?
    Ray is ((xyz origin), (xyz unit direction))
    Voxels are ((xyz),(xyz)...) and assumed 1x1x1"""
    voxels = [v for v in voxel_corners if _intersect(ray,(v,np.add(v,1)))]
    closest_ind = np.argmin([dist_w_offset(ray[0],v,0.5) for v in voxels])
    return voxels[closest_ind]


""" Spiral creation and remembering """

class UnitSpiral:
    """Why a class for this? Because we don't want to recalculate the 
    spiral order each time a spiral gets used, but we want to be able 
    to declare arbitrary numbers of points in the spiral.
    """
    def __init__(self, num_of_pts):
        """Create and remember a spiral with a given num_of_pts"""
        self._n = num_of_pts
        self._xyz, self._rpt = self._fib_sphere(num_of_pts)

    @property
    def n(self):
        """Number of points"""
        return self._n

    @property
    def xyz(self):
        """Cartesian coordinates"""
        return self._xyz

    @property
    def rpt(self):
        """Spherical coordinates (radius, elevation, azimuth)"""
        return self._rpt

    @staticmethod
    def _fib_sphere(n_samples):
        """Sample n points across the surface of a sphere
        
        Points are sampled in cylindrical coords and converted to 
        spherical and cartesian
        """
        s = np.arange(n_samples)
        # Cylindrical theta
        th_0 = 0
        d_th = np.pi * (3 - np.sqrt(5))
        th = np.mod(th_0 + s*d_th, 2*np.pi)
        # Cylindrical z
        d_z = 2/n_samples
        z_0 = d_z/2
        z = (z_0 + s*d_z) - 1
        # Cylindrical r
        r = np.sqrt(1-np.power(z,2))
        # Cartesian x,y. z remains same
        x = np.cos(th) * r
        y = np.sin(th) * r
        xyz = np.stack((x,y,z), 1)
        # Order points using TSP
        distmat = scipy.spatial.distance.squareform(
            scipy.spatial.distance.pdist(xyz))
        first, last = np.argmax(xyz[:,2]), np.argmin(xyz[:,2])
        inds_sorted = greedy.solve_tsp(distmat, 10, endpoints=(first,last))
        xyz = xyz[inds_sorted]
        # Spherical
        rpt = cart_to_sphere(xyz)
        return xyz, rpt





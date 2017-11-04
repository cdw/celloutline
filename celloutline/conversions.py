# encoding: utf-8
""" Raw representation conversions, one type of encoding to the next

We support the following transitions:
 - binary -> (mesh, spiral, spread)
 - spiral -> (point-cloud, mesh)
 - spread -> (binary, mesh)

Author: CDW
"""

# Standard or installed
import numpy as np
import trimesh  #handles mesh creation and interface to OpenSCAD
import skimage.measure
import skimage.segmentation
import scipy.ndimage
import scipy.spatial
# Local
from . import geom


def binary_to_trimesh(binary, step=1):
    """ Convert binary voxel image to a mesh

    Marching cubes meshed output from
    http://scikit-image.org/docs/dev/api/skimage.measure.html#marching-cubes-lewiner

    Parameters
    ----------
    binary: i-by-j-by-k array
        array to mesh, assumed to be binary segmentation
    step: int (default 1)
        number of voxels to step over when generating mesh, larger is courser 
    """
    verts, faces, _, _ = skimage.measure.marching_cubes(binary, step_size=step)
    mesh = trimesh.Trimesh(verts, faces)
    return mesh


def binary_to_spiral(binary, unitspiral=None, num_pts=None):
    """Convert a binary voxel image into a cell spiral

    Parameters
    ---------
    binary: i-by-j-by-k array
        binary voxel segmentation
    unitspiral: None or UnitSpiral object
        unitspiral class to provide rays for intersections
    num_pts: None or int (default 500)
        number of points to use if we need to create our own unitspiral

    Returns
    -------
    spiral_dict: dictionary
        contains the following
        unitspiral: UnitSpiral object
            unitspiral object passed into/created by this function
        radii: num_pts-by-3
            distance to each output point from the origin
        origin: 1-by-3 array
            xyz origin of the spiral
    """
    if num_pts is None:
        num_pts = 500
    if unitspiral is None:
        unitspiral = geom.UnitSpiral(num_pts)
    # Find shell to intersect with
    shell = skimage.segmentation.find_boundaries(
        binary, connectivity=1, mode='outer')
    surf_coords = list(np.array(shell.nonzero()).T)
    surf_voxel_boxes = [(v, np.add(1,v)) for v in surf_coords]
    origin = np.divide(shell.shape,2).astype(int)
    # Find intersections
    xyz = []
    for ray in unitspiral.xyz:
        xyz.append(geom.nearest_intersecting((origin, ray), surf_voxel_boxes))
    radii = [geom.dist_w_offset(origin, pt, 0.5) for pt in xyz]
    spiral_dict = {'unitspiral':unitspiral, 
                   'radii':radii, 
                   'origin':origin}
    return spiral_dict


def binary_to_spread(binary):
    """Convert a binary voxel image into a levelset
    
    The resulting mapping positive inside cell and negative outside

    Parameters
    ----------
    binary: i-by-j-by-k array
        binary voxel segmentation

    Returns
    -------
    spread: i-by-j-by-k array
        a distance mapping where each voxel contains the Euclidean 
        distance to the shell of the segmentation, with a positive sign 
        inside the shell and a negative sign outside the shell
    """
    # The interior is dist to shell
    pos_interior = scipy.ndimage.distance_transform_edt(binary)
    # The exterior is -1 * dist to shell
    neg_exterior = -scipy.ndimage.distance_transform_edt(1-binary)
    # Combine to the levelset
    spread = pos_interior + neg_exterior
    return spread


def spiral_to_point_cloud(radii, unitspiral, origin):
    """Convert a spiral trace back to a xyz point cloud
    
    Parameters
    ----------
    radii: 1-by-n array
        list of distances from the origin to the first shell 
        intersection for a given spiral
    unitspiral: UnitSpiral object
        contains unit rays (in same order as radii) with angles
    origin: 3-by-1 array
        xyz offset of the center of the segmentation

    Returns
    -------
    xyz: n-by-3 array
        locations of each point for which we had a radius
    """
    # Like the original unitspiral, but with new radii
    rpt = np.hstack((np.expand_dims(radii, -1), unitspiral.rpt[:,1:]))
    # Convert the whole thing to cartesian and offset by the origin
    xyz = geom.sphere_to_cart(rpt)
    xyz += origin
    return xyz


def spiral_to_trimesh(radii, unitspiral, origin):
    """Spiral to mesh via a triangulation of the unitspiral's rays"""
    # Convert radii to xyz points
    point_cloud = spiral_to_point_cloud(radii, unitspiral, origin)
    # How would we mesh the original unit-circle spiral?
    faces = scipy.spatial.ConvexHull(unitspiral.xyz).simplices
    # OK, well do that again
    mesh = trimesh.Trimesh(point_cloud, faces)
    # Flip the faces to be all pointing outwards
    mesh.fix_normals()
    return mesh


def spread_to_binary(spread, cutoff):
    """Convert levelset to binary via inequality"""
    return spread>cutoff


def spread_to_trimesh(spread, cutoff):
    """Convert levelset to mesh by way of binary"""
    return binary_to_trimesh(spread_to_binary(spread, cutoff))


# encoding: utf-8
""" Raw representation conversions, one type of encoding to the next

We support the following transitions:
 - binary -> (mesh, spiral, spread)
 - spiral -> (point-cloud, mesh)
 - spread -> (binary, mesh)
 - mesh -> binary

Author: CDW
"""

# Standard or installed
import numpy as np
import trimesh  # handles mesh creation and interface to OpenSCAD
from PIL import Image, ImageDraw # for conversion of mesh to voxel
import skimage.measure
import skimage.segmentation
import scipy.ndimage
import scipy.spatial
# Local
from .utils import geom


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
    mesh.fix_normals()
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
    surf_voxel_boxes = [(v, np.add(1, v)) for v in surf_coords]
    origin = np.divide(shell.shape, 2).astype(int)
    # Find intersections
    xyz = []
    for ray in unitspiral.xyz:
        xyz.append(geom.nearest_intersecting((origin, ray), surf_voxel_boxes))
    radii = [geom.dist_w_offset(origin, pt, 0.5) for pt in xyz]
    spiral_dict = {'unitspiral': unitspiral,
                   'radii': radii,
                   'origin': origin}
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
    rpt = np.hstack((np.expand_dims(radii, -1), unitspiral.rpt[:, 1:]))
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
    return spread > cutoff


def spread_to_trimesh(spread, cutoff):
    """Convert levelset to mesh by way of binary"""
    return binary_to_trimesh(spread_to_binary(spread, cutoff))


def _draw_Path2D(path, shape):
    """Draw a 2d path using Pillow, a support for trimesh_to_binary.

    We are careful to carve out interior voids.
    If the passed path is actually None, return an image containing 0s.

    Parameters
    ----------
    path: trimesh Path2D
        Paths through a plane, generated by intersection through mesh. 
    shape: 3 tuple
        XYZ dimension of output volume

    Returns
    -------
    img: np.array
        Boolean volume
    """
    img = Image.new('1', shape)
    draw = lambda coords, val: ImageDraw.Draw(img).polygon(coords, val)
    # Return zeros for empty paths
    if path is None:
        return np.array(img)
    for poly in path.polygons_full:
        draw(poly.exterior.coords, True)
        for interior in poly.interiors:
            draw(interior.coords, False)
    return np.array(img)


def _take_2d_slice(mesh, z):
    """Take 2D slice of a mesh as a Path2D, support to trimesh_to_binary

    If no intersection occurs at that height, then return None.

    Parameters
    ----------
    mesh: trimesh
        Mesh representing an enclosed volume
    z: float
        z height at which to take section

    Returns
    -------
    section2d: trimesh Path2D
        series of polygons representing the slice
    """
    # Take 3d section
    normal = [0,0,1]
    origin = np.hstack((mesh.centroid[:2], z))
    section3d = mesh.section(normal, origin)
    if section3d is None:
        return None
    # Convert to 2d section
    trans = np.identity(4)
    trans[2,3] = z
    section2d, trans = section3d.to_planar(trans)
    return section2d


def trimesh_to_binary(mesh, shape):
    """Take mesh and original shape, convert mesh to slices"""
    # Where to sample
    z_locs = np.arange(shape[2])+.5 
    # Take sections
    sections = [_take_2d_slice(mesh, z) for z in z_locs]
    # Convert to a binary z-stack 
    z_stack = [_draw_Path2D(path, shape[:2]) for path in sections]
    z_stack = np.swapaxes(np.dstack(z_stack), 0, 1)
    return z_stack

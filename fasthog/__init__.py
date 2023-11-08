"""Reasonably fast implementation of Histogram of Oriented Gradients."""
from ctypes import (CDLL as _CDLL, POINTER as _POINTER, c_double as _c_double,
                    c_void_p as _c_void_p, c_int as _c_int, c_bool as _c_bool)
from ctypes.util import find_library as _find_library
import numpy as _np

_hog_path = _find_library('fasthog')

if not _hog_path:
    import os as _os
    from sys import platform as _platform
    _libroot = _os.path.sep.join(_os.path.realpath(__file__).split(_os.path.sep)[:-5])
    if _platform == 'linux':
        _extension = ".so"
    elif _platform == 'darwin':
        _extension = ".dylib"
    elif _platform == 'win32':
        _extension = '.dll'
    else:
        raise RuntimeError("Invalid operating platform found.")
    _lib = _os.path.join(_libroot, "lib", "libfasthog" + _extension)
    _lib64 = _os.path.join(_libroot, "lib64", "libfasthog" + _extension)
    if _os.path.exists(_lib):
        _hog_path = _lib
    elif _os.path.exists(_lib64):
        _hog_path = _lib64

if not _hog_path:
    raise OSError("Unable to find 'libfasthog'. Add path to its containing directory to your LD_LIBRARY_PATH variable.")

_libhog = _CDLL(_hog_path)

_hog = _libhog.fasthog_hog
_hog.restype = _c_void_p
_hog.argtypes = [
    _POINTER(_c_double),
    _c_int,
    _c_int,
    _c_int,
    _c_int,
    _c_int,
    _c_int,
    _c_int,
    _c_bool,
    _c_int,
    _POINTER(_c_double)
]

_hog_from_gradient = _libhog.fasthog_hog_from_gradient
_hog_from_gradient.restype = _c_void_p
_hog_from_gradient.argtypes = [
    _POINTER(_c_double),
    _POINTER(_c_double),
    _c_int,
    _c_int,
    _c_int,
    _c_int,
    _c_int,
    _c_int,
    _c_int,
    _c_bool,
    _c_int,
    _POINTER(_c_double)
]

_norms = {
    'none': 0,
    'l1': 1,
    'l1_sqrt': 2,
    'l2': 3,
    'l2-hys': 4,
}


def hog(img, cell_size=(8, 8), cells_per_block=(2, 2), n_bins=9, signed=True, norm_type='L2-Hys'):
    """
    Signed histogram of oriented gradients for a given image.

    https://en.wikipedia.org/wiki/Histogram_of_oriented_gradients

    Parameters
    ----------
    img: 2D ndarray dtype=np.float64
         Greyscale image to generate HOG descriptor.
    cell_size: (int, int)
         Cell size in pixels (n_pixels_per_cell_x, n_pixels_per_cell_y). Default: (8, 8).
    cells_per_block: (int, int)
         Number of cells for the block normalization. Default: (2, 2)
    n_bins: int
         Number of bins to histogram orientations. Default: 9.
    signed: bool
         Use signed gradient (0-360) rather than unsigned (0-180). Default: True
    norm_type: str
        'None'
           Don't normalize blocks
        'L1'
           Normalization using L1-norm.
        'L1-sqrt'
           Normalization using L1-norm, followed by square root.
        'L2'
           Normalization using L2-norm.
        'L2-Hys'
           Normalization using L2-norm, followed by limiting the
           maximum values to 0.2 (`Hys` stands for `hysteresis`) and
           renormalization using L2-norm. (default)

    Returns
    -------
    np.ndarray((n_blocks_x, n_blocks_y, n_bins), dtype=np.float64)
         Hog descriptor for image. Right/bottom sides of image are truncated if image is not an integer
         multiple of the cell size. I.e.
         n_cells_x = img.shape[1] // cell_size[1]
         n_cells_y = img.shape[0] // cell_size[0]
         n_blocks_x = (n_cells_x - cells_per_block[0]) + 1
         n_blocks_y = (n_cells_y - cells_per_block[1]) + 1
    """
    norm_type = norm_type.lower()
    norm_index = _norms.get(norm_type, None)
    if norm_index is None:
        print(f"Invalid norm type for hog: {norm_type}")
        return None
    n_cells_x = img.shape[1] // cell_size[1]
    n_cells_y = img.shape[0] // cell_size[0]
    n_blocks_y = (n_cells_y - cells_per_block[1]) + 1
    n_blocks_x = (n_cells_x - cells_per_block[0]) + 1
    res = _np.empty((n_blocks_y, n_blocks_x, n_bins))
    _hog(img.ctypes.data_as(_POINTER(_c_double)),
         img.shape[1],
         img.shape[0],
         cell_size[0],
         cell_size[1],
         cells_per_block[0],
         cells_per_block[1],
         n_bins,
         signed,
         norm_index,
         res.ctypes.data_as(_POINTER(_c_double)))
    return res


def hog_from_gradient(gx, gy, cell_size=(8, 8), cells_per_block=(2, 2), n_bins=9, signed=True, norm_type='L2-Hys'):
    """
    Signed histogram of oriented gradients given the gradient of the image.

    https://en.wikipedia.org/wiki/Histogram_of_oriented_gradients

    Parameters
    ----------
    gx: 2D ndarray dtype=np.float64
         gradient of image w.r.t. x
    gy: 2D ndarray dtype=np.float64
         gradient of image w.r.t. y
    cell_size: (int, int)
         Cell size in pixels (n_pixels_per_cell_x, n_pixels_per_cell_y). Default: (8, 8).
    cells_per_block: (int, int)
         Number of cells for the block normalization. Default: (2, 2)
    n_bins: int
         Number of bins to histogram orientations. Default: 9.
    signed: bool
         Use signed gradient (0-360) rather than unsigned (0-180). Default: True
    norm_type: str
        'None'
           Don't normalize blocks
        'L1'
           Normalization using L1-norm.
        'L1-sqrt'
           Normalization using L1-norm, followed by square root.
        'L2'
           Normalization using L2-norm.
        'L2-Hys'
           Normalization using L2-norm, followed by limiting the
           maximum values to 0.2 (`Hys` stands for `hysteresis`) and
           renormalization using L2-norm. (default)


    Returns
    -------
    np.ndarray((n_blocks_x, n_blocks_y, n_bins), dtype=np.float64)
         Hog descriptor for image. Right/bottom sides of image are truncated if image is not an integer
         multiple of the cell size. I.e.
         n_cells_x = img.shape[1] // cell_size[1]
         n_cells_y = img.shape[0] // cell_size[0]
         n_blocks_x = (n_cells_x - cells_per_block[0]) + 1
         n_blocks_y = (n_cells_y - cells_per_block[1]) + 1
    """
    n_cells_x = gx.shape[1] // cell_size[1]
    n_cells_y = gx.shape[0] // cell_size[0]
    n_blocks_y = (n_cells_y - cells_per_block[1]) + 1
    n_blocks_x = (n_cells_x - cells_per_block[0]) + 1
    norm_type = norm_type.lower()
    norm_index = _norms.get(norm_type, None)
    if norm_index is None:
        print(f"Invalid norm type for hog: {norm_type}")
        return None

    res = _np.empty((n_blocks_y, n_blocks_x, n_bins))
    _hog_from_gradient(gx.ctypes.data_as(_POINTER(_c_double)),
                       gy.ctypes.data_as(_POINTER(_c_double)),
                       gx.shape[1],
                       gx.shape[0],
                       cell_size[0],
                       cell_size[1],
                       cells_per_block[0],
                       cells_per_block[1],
                       n_bins,
                       signed,
                       norm_index,
                       res.ctypes.data_as(_POINTER(_c_double)))
    return res

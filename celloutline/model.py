# encoding: utf-8
""" Learned dimensional reductions
Author: CDW
"""

# Standard or installed
import pickle
import numpy as np
from sklearn.decomposition import IncrementalPCA
# Local
from . import representation

__all__ = ["PCA", "ICA", "Fourier", "Autoencoder"]


class Model:
    """Model is a base class that provides some "don't repeat yourself"
    for models which inherit from it.
    """

    def __init__(self, existing_fn=None):
        pass

    def load(self, fn):
        """Load an existing model"""
        if not fn.endswith('.pkl'):
            fn += '.pkl'
        self.fn = fn
        with open(fn, 'rb') as modelfile:
            self.model = pickle.load(modelfile)

    def save(self, fn=None):
        """Save model to file, existing or new (if fn is passed)"""
        if fn is None:
            fn = self.fn
        elif not fn.endswith('.pkl'):
            fn += '.pkl'
        with open(fn, 'wb') as modelfile:
            pickle.dump(self.model, modelfile)

    def set_params(self, params=None):
        """Set parameters from passed or local"""
        if params is not None:
            self.model.set_params(**params)
        else:
            self.model.set_params(**self.params)


class PCA(Model):
    """Given a set of input vectors, find their principle components"""

    def __init__(self, fn=None, n_comp=None, batch_size=None):
        self.model = IncrementalPCA()
        self.fn = fn
        self.params = {"n_components": n_comp, "batch_size": batch_size}
        self.set_params()

    def load(self, fn):
        """Set parameters after loading from filename"""
        super().load(fn)
        self.params = self.model.get_params()
        return

    def fit(self, reps):
        """Fit a list of representations"""
        X = [r.to_vector() for r in reps]
        self.model.fit(X)

    def err(self, to_transform, to_check_against):
        """Mesh error between reconstructed to_transform representation and
        mesh conversion of to_check_against
        """
        vec = to_transform.to_vector()
        vec_trans = self.model.transform(vec)
        vec_recon = self.model.inverse_transform(vec_trans)
        transformed = to_transform.from_vector(vec_recon)
        mesh1 = transformed.mesh()
        mesh2 = to_check_against.mesh()
        error = representation.mesh_error(mesh1, mesh2)
        return error


class ICA:
    def __init__(self, n_components):
        """Reduce/replicate an image using independent components

        Parameters
        ----------
        n_components: int
            number of independent components to retain
        """
        self.n_components = n_components

    def fit(self):
        pass

    def error(self, original, recon):
        pass


class Fourier:
    def __init__(self, reduced_dimensions):
        """Reduce/replicate an image using its Fourier bases

        Parameters
        ----------
        reduced_dimensions: int
            number of dimensions to keep when encoding,
            is coerced into being an even cube if needed.
            Even cubes are (2*n)**3, e.g. 8, 64, 216, 512,
            1000, 1728, 2744, 4096, 5832, 8000, ...
        """
        self.reduced_dimensions = reduced_dimensions

    @property
    def reduced_dimensions(self):
        """How many dimensions to keep, see __init__"""
        return self._reduced_dimensions

    @reduced_dimensions.setter
    def reduced_dimensions(self, dims):
        """Enforce even cube, see __init__ for documentation"""
        # Coerce reduced dimensions into even cube form
        to_keep = round(np.divide(np.power(dims, 1 / 3), 2))
        cubed_dims = np.power(np.multiply(to_keep, 2), 3).astype(int)
        self._reduced_dimensions = cubed_dims
        # Give a warning if the coercion changed the number of dims
        if cubed_dims != dims:
            import warnings
            warnings.warn("Passed dimensions, %i, not an even cube. \
                          Converted to %i" % (reduced_dimensions, cubed_dims))
        return

    @property
    def _mask(self):
        """Derive mask from original shape"""
        center = np.divide(self.original_shape, 2).astype(int)
        to_keep = int(
            round(np.divide(np.power(self.reduced_dimensions, 1 / 3), 2)))
        mask = np.zeros(self.original_shape)
        mask[center[0] - to_keep:center[0] + to_keep, center[1] - to_keep:
             center[1] + to_keep, center[2] - to_keep:center[2] + to_keep] = 1
        return mask

    @property
    def original_shape(self):
        """Original image or input shape"""
        return self._original_shape

    @original_shape.setter
    def original_shape(self, input):
        """Set the original shape, input must have .shape property"""
        self._original_shape = input.shape

    def encoded_to_flat(self, encoded):
        """Flatten the encoded representation"""
        assert encoded.shape == self.original_shape, "Shape mismatch, reset"
        selected_freqs = np.fft.fftshift(encoded)[self._mask.astype(bool)]
        return selected_freqs.flatten()

    def flat_to_encoded(self, flat):
        """Return flat representation to encoded dimensions for decoding"""
        cube_side = int(round(np.power(self.reduced_dimensions, 1 / 3)))
        dimensions = [cube_side for dimension in self.original_shape]
        selected_freqs = flat.reshape(dimensions)
        to_pad = [((dimension - cube_side) // 2, (dimension - cube_side) // 2)
                  for dimension in self.original_shape]
        encoded = np.pad(selected_freqs, to_pad, 'constant')
        return encoded

    def encode(self, input):
        """Encode passed image

        Parameters
        ----------
        input: numpy array
            3D numpy array
        """
        # Take n dimensional fft of input
        fftcell = np.fft.fftn(input)
        # Create mask to keep just n frequencies
        self.original_shape = input
        mask = self._mask
        # Mask off only ~n frequencies
        shiftedfftcell = np.fft.fftshift(
            fftcell)  # shift max freq around center
        fftless = shiftedfftcell * mask  # mask
        return fftless

    def decode(self, encoded):
        """Decode passed encoded form"""
        ifftless = np.fft.ifftn(np.fft.ifftshift(encoded))
        return ifftless


class Autoencoder:
    def __init__(self, train, test):
        pass

    def fit(self):
        pass

    def err(self, original, recon):
        pass

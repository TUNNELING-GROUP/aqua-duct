from distutils.core import setup
from Cython.Build import cythonize


setup(
    ext_modules = cythonize('aquarium/*.py')
)
setup(
    ext_modules = cythonize('aquarium/*/*.py')
)
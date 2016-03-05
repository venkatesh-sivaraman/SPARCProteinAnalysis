from distutils.core import setup, Extension
from Cython.Build import cythonize

setup (name = "PythonProteins",
	   version = "1.0",
	   description = "Protein folding simulation tool",
	   ext_modules = cythonize("proteinmath.pyx"))
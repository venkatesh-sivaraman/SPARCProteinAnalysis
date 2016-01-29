from distutils.core import setup, Extension

module1 = Extension('fastint', ['_fastint.c'])

setup (name = "stringint",
	   version = "1.0",
	   description = "Fast conversion of string to int.",
	   ext_modules = [module1])
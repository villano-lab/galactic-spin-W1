# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
# import os
# import sys
# sys.path.insert(0, os.path.abspath('.'))
import sphinx_rtd_theme
import subprocess, os, sys
sys.path.append('../binder/python')
#real quick, make a link if it doesn't exist
os.chdir('../')
try:
    os.symlink('../binder/data','data')
except FileExistsError:
    pass
os.chdir('source')

autodoc_mock_imports = ['matplotlib','IPython','astroquery','astropy','ipywidgets']

# Check if we're running on Read the Docs' servers
read_the_docs_build = os.environ.get('READTHEDOCS', None) == 'True'

breathe_projects = {
    "galactic-spin-W1":"xml",
}

if read_the_docs_build:
    input_dir = '../../binder/python'
    output_dir = 'dox_build'
    sys.path.append('../../binder/python')
    #os.symlink('../../binder/data','data')
    #os.symlink('../../binder/python','python')


# -- Project information -----------------------------------------------------

project = 'galactic-spin-W1'
copyright = '2022, A.N. Villano, J. Bergfalk, K. Harris, F. Vititoe, R. Hatami, J. Johnston'
author = 'A.N. Villano, J. Bergfalk, K. Harris, F. Vititoe, R. Hatami, J. Johnston'

# The full version, including alpha/beta/rc tags
release = '2.0.2'

# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = ['sphinx_rtd_theme','myst_parser','sphinx.ext.autodoc','sphinx.ext.autosummary','autoapi.extension']
#autodoc_default_flags = ['members']
autosummary_generate = True
autoapi_dirs = ['../../binder/python']
autoapi_template_dir = 'autoapi_templates'
autoapi_keep_files = True #helps debugging

# Breathe Configuration
#breathe_default_project = "galactic-spin-W1"

# Add any paths that contain templates here, relative to this directory.
templates_path = ['templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ['build','templates','autoapi_templates']
autoapi_ignore = ['*.cpython-??.pyc','*-checkpoint.py*']


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
#html_theme = 'alabaster'
html_theme = 'sphinx_rtd_theme'

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = []
==================
1. Getting started
==================

Binder
======

.. image:: https://mybinder.org/badge_logo.svg
   :target: https://mybinder.org/v2/gh/villano-lab/galactic-spin-W1/master?labpath=binder

Binder is a website that hosts jupyter notebooks and allows users to run them on the web without any download or installation required
from the user. If you only wish to run the notebooks, you can click on the binder button above to run them from your browser.

Installation
============

Installation of this workshop is not needed. 
Binder Project provides an easy access to Jupyter notebooks and their environments without the hassle of installation. 
Upon clicking on the provided Binder URLs, the modules of the workshop open in your browser. 

If you wish to work offline, you can run these notebooks in Jupyter Lab or Jupyter Notebooks.

Requirements
------------

#. **Jupyter Lab** or **Jupyter Notebooks**, which are both included with `Anaconda Python <https://www.anaconda.com/>`_.
#. The ability to install python modules, either through **pip**, **conda**, or `Anaconda Python <https://www.anaconda.com/>`_.
#. The following python modules available through conda, as listed in `binder/environment.yml <https://github.com/villano-lab/galactic-spin-W1/tree/master/binder/environment.yml>`_:
   
   * python 3.7.9
   
   * numpy
   
   * scipy
   
   * matplotlib 3.5.1
   
   * mpld3
   
   * sympy
   
   * ipympl
   
   * ipylab
   
   * ipywidgets 8.0.4
   
   * ffmpeg
   
   * ipykernel
   
   * jupyterlab 3.4.3
   
   * pip

#. The following python modules available through pip, as listed in `binder/environment.yml <https://github.com/villano-lab/galactic-spin-W1/tree/master/binder/environment.yml>`_:
   
   * lmfit
   
   * astroquery
   
   * dataPython

Instructions
------------

These instructions will assume that you will use Anaconda. 
For experienced users with other setups, 
please install the environment found in `binder/environment.yml <https://github.com/villano-lab/galactic-spin-W1/tree/master/binder/environment.yml>`_ to run our notebooks.

0. Install `Anaconda Python <https://www.anaconda.com/>`_.
1. `Download <https://github.com/villano-lab/galactic-spin-W1/archive/refs/heads/master.zip>`_ and unzip the repository somewhere inside your user directory, such as the "Documents" folder.
2. In Anaconda Navigator, go to "Environments," then select "Import."
3. Import from "Local drive" and select the file `binder/environment.yml <https://github.com/villano-lab/galactic-spin-W1/tree/master/binder/environment.yml>`_ from the repository. Name the environment `galactic-spin`. Anaconda will start installing all other necessary dependencies; please be patient, as this will take some time.
4. Once the environment is finished building, from the "Home" menu of Anaconda Navigator, launch a prompt. The prompt will depend on your system -- for example, in Windows, you can choose "CMD.exe Prompt."
5. In the prompt, run the following command: `ipython kernel install --user --name=galactic-spin`
6. In Anaconda Navigator, launch "JupyterLab" and browse to your local copy of the repository, then enter the `binder` subdirectory. Enjoy!

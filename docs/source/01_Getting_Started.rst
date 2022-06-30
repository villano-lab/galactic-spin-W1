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

1. **Jupyter Lab** or **Jupyter Notebooks**, which are both included with `Anaconda Python <https://www.anaconda.com/>`_.
2. The ability to install python modules, either through **pip**, **conda**, or `Anaconda Python <https://www.anaconda.com/>`_.

Instructions
------------

These instructions will assume that you will use Anaconda. 
For experienced users with other setups, 
please install the environment found in `binder/environment.yml` to run our notebooks.

0. Install `Anaconda Python <https://www.anaconda.com/>`_.
1. `Download <https://github.com/villano-lab/galactic-spin-W1/archive/refs/heads/master.zip>`_ and unzip the repository somewhere inside your user directory, such as the "Documents" folder.
2. In Anaconda Navigator, go to "Environments," then select "Import."
3. Import from "Local drive" and select the file `binder/environment.yml` from the repository. Name the environment `galactic-spin`. Anaconda will start installing all other necessary dependencies; please be patient, as this will take some time.
4. Once the environment is finished building, from the "Home" menu of Anaconda Navigator, launch a prompt. The prompt will depend on your system -- for example, in Windows, you can choose "CMD.exe Prompt."
5. In the prompt, run the following command: `ipython kernel install --user --name=galactic-spin`
6. In Anaconda Navigator, launch "JupyterLab" and browse to your local copy of the repository, then enter the `binder` subdirectory. Enjoy!

Doxygen test
============

.. doxygenfunction:: galdict
   :project: galactic-spin-W1
   :path: xml
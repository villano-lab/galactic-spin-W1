import subprocess
import tempfile
import os
from pytest_notebook.execution import execute_notebook
from pytest_notebook.notebook import load_notebook

n = [0,0,0,0,0,0,0,0,0,0]
n[0] = load_notebook('01_DM_Rotation_Curve_Intro.ipynb')
n[1] = load_notebook('02_Widget_NGC5533_DMonly.ipynb')
n[2] = load_notebook('03_Measured_Data_Plotting.ipynb')
n[3] = load_notebook('04_Plotting_Rotation_Curves.ipynb')
n[4] = load_notebook('05_Widget_NGC5533_All_Components.ipynb')
n[5] = load_notebook('06_Plotting_SPARC_data.ipynb')
n[6] = load_notebook('07_Bonus_Bulge_Rotation_Curve.ipynb')
n[7] = load_notebook('08_Interactive_Fitting.ipynb')
n[8] = load_notebook('09_Widget_SPARC_Galaxies.ipynb')
n[9] = load_notebook('10_Bonus_Black_Holes_as_DM.ipynb')

def test():
    print('Testing Jupyter notebooks...')
    for notebook in n:
        execute_notebook(notebook,with_coverage=True)

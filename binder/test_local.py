import subprocess
import tempfile
import os

def _exec_notebook(path):
    with tempfile.NamedTemporaryFile(suffix=".ipynb") as fout:
        args = ["jupyter", "nbconvert", "--to", "notebook", "--execute",
                "--ExecutePreprocessor.timeout=10000",
                "--ExecutePreprocessor.kernel_name=python3",
                "--output", fout.name, path]
        subprocess.check_call(args)


def test():
    print('Testing Jupyter notebooks...')
    _exec_notebook('01_DM_Rotation_Curve_Intro.ipynb')
    _exec_notebook('02_Widget_NGC5533_DMonly.ipynb')
    _exec_notebook('03_Measured_Data_Plotting.ipynb')
    _exec_notebook('04_Plotting_Rotation_Curves.ipynb')
    _exec_notebook('05_Widget_NGC5533_All_Components.ipynb')
    _exec_notebook('06_Plotting_SPARC_data.ipynb')
    _exec_notebook('07_Bonus_Bulge_Rotation_Curve.ipynb')
    _exec_notebook('08_Interactive_Fitting.ipynb')
    _exec_notebook('09_Widget_SPARC_Galaxies.ipynb')
    _exec_notebook('10_Bonus_Black_Holes_as_DM.ipynb')

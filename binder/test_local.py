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
    _exec_notebook('1_DM_Rotation_Curve_Intro.ipynb')
    _exec_notebook('2_Measured_Data_Plotting.ipynb')
    _exec_notebook('3_Plotting_Rotation_Curves.ipynb')
    _exec_notebook('4_Interactive_Fitting.ipynb')
    _exec_notebook('5_Plotting_SPARC_data.ipynb')
    _exec_notebook('Bonus1_Bulge_Decomposition.ipynb')
    _exec_notebook('Bonus2_Black_Holes_as_DM.ipynb')
    _exec_notebook('Widget1_NGC5533_DM_only.ipynb')
    _exec_notebook('Widget2_NGC5533_All_Components.ipynb')
    _exec_notebook('Widget3_SPARC_Galaxies.ipynb')

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
    _exec_notebook('Bulge_Activity.ipynb')
    _exec_notebook('Interactive_Fitting.ipynb')
    _exec_notebook('Interactive_Measured_Data_Plotting.ipynb')
    _exec_notebook('Interactive_Rotation_Curve_Plotting.ipynb')
    _exec_notebook('NGC5533_DM_only_widget.ipynb')
    #_exec_notebook('RC_Sliders-Multiple_Galaxies.ipynb') 
    _exec_notebook('Rotation_Curve_Animation.ipynb')
    _exec_notebook('SPARC_activity_newtonian_mass_models.ipynb')

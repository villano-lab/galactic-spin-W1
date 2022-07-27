"""A module for handling the downloading of SPARC data in notebook `09_Widget_SPARC_Galaxies.ipynb <https://github.com/villano-lab/galactic-spin-W1/blob/master/binder/09_Widget_SPARC_Galaxies.ipynb>`_.
"""

###########################
### Download SPARC data ###
###########################

import numpy as np
import requests
import IPython

def downloadsparc():
    """
    Downloads rotation curve data for Newtonian Mass Models (Rotmod_LTG.zip) from the SPARC website: http://astroweb.cwru.edu/SPARC/.

    Parameters:
        None 

    Returns:
        None. Files are downloaded to local computer.

    Example:
        >>> downloadsparc()
        Download URL: http://astroweb.cwru.edu/SPARC/Rotmod_LTG.zip.
        Downloaded zip file successfully.
    """
    
    out.clear_output(True)
    url = 'http://astroweb.cwru.edu/SPARC/Rotmod_LTG.zip'         # Define url where the original zip file is located
    myfile = requests.get(url, allow_redirects=True)              # Download file
    open('data/sparc/Rotmod_LTG.zip', 'wb').write(myfile.content) # Define location and filename where to download file
    print("Download URL: {}.".format(url))
    print("Downloaded zip file successfully.")
    

##########################
####### Unzip file #######
##########################

import zipfile
import pathlib
import os

def unzipfiles():
    """
    Extracts all files of the zipped rotation curve data (Rotmod_LTG.zip), downloaded from the SPARC website: http://astroweb.cwru.edu/SPARC/.

    Parameters:
        None 

    Returns:
        None. Unzips files on local computer to the 'data/sparc/' folder.

    Example:
        >>> downloadsparc()
        Download URL: http://astroweb.cwru.edu/SPARC/Rotmod_LTG.zip.
        Downloaded zip file successfully.
        >>> unzipfiles()
        Data files unzipped successfully.
        SPARC data file location: '/mnt/c/Users/harrkath/Documents/GitHub/galactic-spin-W1/binder/data/sparc/'.
    """
    
    with zipfile.ZipFile('data/sparc/Rotmod_LTG.zip', 'r') as zip_ref:    # Define zip file location and name
        location = 'data/sparc/'                                          # Define location where to unzip files
        zip_ref.extractall(location)                                      # Extract all
    print("Data files unzipped successfully.")
    path = str(pathlib.Path(__file__).parent.resolve())                   # Find the path of the file 
    path = path.replace("python","")                                      # Data files are not in the python folder
    path = path.replace(os.sep,"/")                                       # Data files are not in the python folder
    print("SPARC data file location: '{}'.".format(path + location))
    

#################################
########### Buttons #############
#################################

from ipywidgets import Button, HBox, VBox, Output
from IPython.display import Javascript, display

buttonYES = Button(                   # YES button
    description="Yes",
    button_style='success', # 'success', 'info', 'warning', 'danger' or ''
    icon='check')
"""Button, WITHOUT display box, for downloading and unzipping SPARC data.

:type: ipywidgets.widgets.widget_button.Button
"""

out = Output()
"""Displaying the button.

:type: ipywidgets.widgets.Output
"""

def on_button_clicked_YES(_):         # When clicked on the 'YES' button, download and unzip SPARC data, run all cells below
    """
    Button labeled as 'YES' to download and unzip SPARC data.

    Parameters:
        None 

    Returns:
        None. Downloads files to local computer and unzips them.
    """
    
    with out:
        out.clear_output()
        downloadsparc()
        unzipfiles()
        #display(Javascript('IPython.notebook.execute_cell_range(IPython.notebook.get_selected_index()+1, IPython.notebook.ncells())'))   #not working

buttonYES.on_click(on_button_clicked_YES)

displaybuttons = VBox([HBox([buttonYES]),out])
"""Visible button for downloading and unzipping SPARC data.

Version of :func:`buttonYES <downloadSPARCdata.buttonYES>` that can be displayed without issue.

:type: ipywidgets.widgets.widget_box.VBox
"""


#####################################
### Choose a galaxy from filelist ###
#####################################

from pprint import pprint
from ipywidgets import Dropdown, Widget, Box, Layout, Label
import ipywidgets as w

# Make a list of galaxies from sparc data filenames
#filelist_raw = os.listdir('./data/sparc/')                   
#filelist = [f for f in os.listdir('./data/sparc/') if "rotmod.dat" in f]    # Only look at filenames that has 'rotmod.dat' in the name
options = []
"""A list of galaxies to choose from.

:type: list

.. seealso:: This list is generated based on the files present in the :code:`./data/sparc` directory, which is populated by the :func:`downloadsparc <downloadSPARCdata.downloadsparc>` and :func:`unzipfiles <downloadSPARCdata.unzipfiles` functions.
"""

for f in [f for f in os.listdir('./data/sparc/') if "rotmod.dat" in f]:
    #galaxy = f.replace("_rotmod.dat","")                     # Remove "_rotmod.dat" from filename
    options.append(f.replace("_rotmod.dat",""))

# Create a dropdown menu of all galaxies
form_item_layout = Layout(
    display='flex',
    flex_flow='row',
    justify_content='space-between'
)
"""A layout to be used by :func:`galaxyoptions <downloadSPARCdata.galaxyoptions>` for the :func:`dropdown <downloadSPARCdata.dropdownmenu>`.

:type: ipywidgets.widgets.widget_layout.Layout
"""

dropdownmenu = Dropdown(options=options, value='NGC5005')
"""A dropdown menu for selecting a galaxy.

:type: ipywidgets.widgets.widget_selection.Dropdown
"""

galaxyoptions = Box([Box([Label(value='Galaxy: '),dropdownmenu], layout=form_item_layout)], layout=Layout(
    display='flex',
    flex_flow='column',
    border='solid 2px',
    align_items='stretch',
    width='50%'
))
"""A displayable box, including :func:`a dropdown menu <downloadSPARCdata.dropdownmenu>`, for selecting a galaxy.

:type: ipywidgets.widgets.widget_box.Box
"""
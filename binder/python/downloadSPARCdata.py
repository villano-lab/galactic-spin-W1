###########################
### Download SPARC data ###
###########################

import numpy as np
import requests
import IPython

def downloadsparc():
    """
    Downloads Rotmod_LTG.zip from the SPARC website: http://astroweb.cwru.edu/SPARC/.

    Parameters:
        None 

    Returns:
        File downloaded on local computer.
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
    Extracts all files of the Rotmod_LTG.zip, dowloaded from the SPARC website: http://astroweb.cwru.edu/SPARC/.

    Parameters:
        None 

    Returns:
        Unzipped files on local computer in the 'data/sparc/' folder.
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
out = Output()

def on_button_clicked_YES(_):         # When clicked on the 'YES' button, download and unzip SPARC data, run all cells below
    """
    'YES' button to download and unzip SPARC data.

    Parameters:
        None 

    Returns:
        None.
    """
    
    with out:
        out.clear_output()
        downloadsparc()
        unzipfiles()
        #display(Javascript('IPython.notebook.execute_cell_range(IPython.notebook.get_selected_index()+1, IPython.notebook.ncells())'))   #not working

buttonYES.on_click(on_button_clicked_YES)

buttons = HBox([buttonYES])
displaybuttons = VBox([buttons,out])


#####################################
### Choose a galaxy from filelist ###
#####################################

from pprint import pprint
from ipywidgets import Dropdown, Widget, Box, Layout, Label
import ipywidgets as w

# Make a list of galaxies from sparc data filenames
filelist_raw = os.listdir('./data/sparc/')                   
filelist = [f for f in filelist_raw if "rotmod.dat" in f]    # Only look at filenames that has 'rotmod.dat' in the name
options = []
for f in filelist:
    galaxy = f.replace("_rotmod.dat","")                     # Remove "_rotmod.dat" from filename
    options.append(galaxy)
options = np.array(options)

# Create a dropdown menu of all galaxies
form_item_layout = Layout(
    display='flex',
    flex_flow='row',
    justify_content='space-between'
)

dropdownmenu = Dropdown(options=options, value='NGC5005')

form_items = [
    Box([Label(value='Galaxy: '),dropdownmenu
         ], layout=form_item_layout)
]

galaxyoptions = Box(form_items, layout=Layout(
    display='flex',
    flex_flow='column',
    border='solid 2px',
    align_items='stretch',
    width='50%'
))
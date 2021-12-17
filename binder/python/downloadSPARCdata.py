###########################
### Download SPARC data ###
###########################

import numpy as np
import requests

def downloadsparc():
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
buttonNO = Button(                    # NO button
    description="No",
    button_style='warning', # 'success', 'info', 'warning', 'danger' or ''
    icon='close')
out = Output()

def on_button_clicked_YES(_):         # When clicked on the 'YES' button, download and unzip SPARC data, run all cells below
    downloadsparc()
    unzipfiles()
    display(Javascript('IPython.notebook.execute_cell_range(IPython.notebook.get_selected_index()+1, IPython.notebook.ncells())'))   

def on_button_clicked_NO(_):          # When clicked on the 'NO' button, run all cells below
    display(Javascript('IPython.notebook.execute_cell_range(IPython.notebook.get_selected_index()+1, IPython.notebook.ncells())'))         

buttonYES.on_click(on_button_clicked_YES)
buttonNO.on_click(on_button_clicked_NO)

buttons = HBox([buttonYES,buttonNO])
displaybuttons = VBox([buttons,out])


#####################################
### Choose a galaxy from filelist ###
#####################################

from pprint import pprint
from ipywidgets import Dropdown, Widget
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
galaxyoptions = Dropdown(                                    
    options=options,
    value='NGC5005',
    description='Galaxy:',
    disabled=False,
)
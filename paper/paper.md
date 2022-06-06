---
title: 'The Data Behind Dark Matter: Exploring Galactic Rotation' 
bibliography: references.bib
tags:
  - Python
  - Jupyter 
  - Binder
  - Workshop 
  - Astronomy
authors:
  - name: A.N. Villano
    affiliation: 1
    orcid: 0000-0002-3893-7259
  - name: Kitty C. Harris
    affiliation: 2
    orcid: 0000-0001-5406-8367
  - name: Judit Bergfalk 
    affiliation: 3
    orcid: 0000-0003-1662-0768 
  - name: Raphael Hatami
    affiliation: 1
  - name: F. Vititoe 
    affiliation: 4
  - name: Julia Johnston 
    affiliation: 3
affiliations:
 - name: Department of Physics, University of Colorado Denver, Denver CO 80217, USA
   index: 1
 - name: Integrated Sciences, University of Colorado Denver, Denver CO 80217, USA
   index: 2
 - name: Astrophysical & Planetary Sciences, University of Colorado Boulder, Boulder, CO 80309, USA 
   index: 3
 - name:  Department of Physics, University of Colorado Boulder, Boulder, CO 80309, USA 
   index: 4
date: 30 March 2022
nocite: '@*'
---

# Summary

By analyzing the rotational velocities of bodies in galaxies, physicists and astronomers have
found that there seems to be something missing in our understanding of these galaxies. One theory
is that there is some invisible matter present that does not interact with light
- that is to say, these galaxies contain dark matter [@1978ApJ-225L-107R]. 

Participants in this workshop will have the opportunity to explore dark matter through scientific
literature-based [@Karukes_2015;@Richards_2015;@Fraternali_2011;@de_Naray_2008] galactic rotation
curves both by using interactive programs and by editing Python code. This will give participants
an understanding of how physicists arrived at the idea of dark matter, showing them the difference
between curve fits with and without dark matter components. Understanding dark matter's
epistemological origins will help participants formulate their own opinions on the dark matter
debate.

## Description of the software or learning module

In this interactive module, we include several modules in the form of Jupyter notebooks [@soton403913]:

| File Name | Short Description |
| --- | --- |
| 01_DM_Rotation_Curve_Intro.ipynb | Animations and rotation curve plots demonstrating three types of rotational motion. |
| 02_Widget_NGC5533_DMonly.ipynb | Interactive introduction to dark matter. | 
| 03_Measured_Data_Plotting.ipynb | Rotation curve plotting of measured velocities to visualize star and gas motions in a galaxy. | 
| 04_Plotting_Rotation_Curves.ipynb | Plotting the rotation curves of galaxy components. | 
| 05_Widget_NGC5533_All_Components.ipynb | Interactive widget to visualize the components of the galaxy NGC 5533. | 
| 06_Plotting_SPARC_Data.ipynb | Plotting the components of galactic rotation curves using the SPARC database of 175 galaxies. | 
| 07_Bonus_Bulge_Rotation_Curve.ipynb| Constructing a rotation curve for the bulge component using empirically-derived parameters. | 
| 08_Interactive_Fitting.ipynb | Interactive curve fitting. | 
| 09_Widget_SPARC_Galaxies.ipynb | Interactive widget to visualize the components of multiple galaxies using the SPARC database. | 
| 10_Bonus_Black_Holes_as_DM.ipynb | Considering tiny black holes as dark matter candidates. | 


# Statement of Need 

Rotation curves are a key empirical artifact through which dark matter can be observed and analyzed [@1978ApJ-225L-107R].
However, a thorough, start-to-finish description of the rotation curve building process is typically not given in scientific publications. Furthermore, software tools used in rotation curve literature are generally difficult for inexperienced users; 
for example, the GIPSY software package is very thorough but does not provide any introduction 
as it is intended for experienced users with a firm grasp on rotation curve components [@Gipsy_1992]. 
Therefore, a rigorous yet accessible learning module is needed to provide an entry point 
for any individual interested in investigating the effect of dark matter in spiral galaxies. 
Our workshop is designed to present a convenient platform for developing basic rotation curves 
focused on introducing newcomers to the concepts necessary for understanding galactic rotation. 
This is acheived by leading users through hands-on computational activities, 
including building and plotting their own rotation curves. 
The primary goal of our project is to present rotation curve development and research in a versatile and approachable format 
for anyone to explore, learn from, and build upon.

# Learning Objectives 

The learning objectives for this educational module are:
    1. Provide a working space where people can connect with current literature and identify as
scientists.  
    2. Educate curious students or other individuals on the basic concepts of rotation curves, as related to the current problems and mysteries regarding dark matter in the universe.
    3. Provide users with accessible activities relating to the basic principles of rotation curve composition. This includes:
        a. facilitating the introduction of rotation curve concepts via open-source code.
        b. interactive programs to provide users with practical and tangible approach of what producing rotation curves involves.
    4. Learn to use the SPARC database to plot rotation curves of many galaxies.
    5. Understand data and models by interacting directly with equations and figures.

Most of the content provided in this module has been presented and taught in previous
workshops/research symposiums (University of Colorado Denver: Data Science Symposium 2021,
Research and Creative Activities Symposium 2020, 2021, and 2022 [REFS]) with feedback collected from
participants. We have chosen the activities for this module that proved most successful in terms
of education and sparking interest in participants. 

# Delivery

The modules are designed to be presented to students in numeric order as part of a workshop, 
skipping those marked as "Bonus" as needed to fit the alloted time. 
Students are encouraged to work together to complete the modules and compare their results to one another. 
While working through a module, the instructor(s) should be available to answer questions and check in on participants' progress, 
but they should leave the bulk of the work to the participants themselves. If any bonus modules are being skipped, 
the instructor(s) may wish to suggest these to participants who find they are completing the content at a faster pace than the rest of the room.

# Story 

This project emerged from several years of literary analysis, trial and error, and
reproducibility studies our team conducted on rotation curve research. A principal impression our
team drew from this journey was the lack of clarity and accessibility in the world of rotation
curve research, including programs and data requisite for composing rotation curves. Examining the
difficulties encountered in our team's research, we developed an educational module that is aimed
to be clear and concise, with curve-composition activities that can be easily reproduced by the
user. From our experience analyzing a number of rotation curve publications in addition to
correspondences with active researchers in the field and with currently-used curve-composing
programs, we have begun to develop our own versions of such programs, making improvements which we
believe will greatly improve accessibility and understanding. We also include content which we
produced for previous workshops and research symposiums that further facilitate the understanding
of producing rotation curves. 

# Acknowledgements

The authors would like to thank Dr. Martin Vogelaar at Kapteyn Astronomical Institute, Dr. Edo
Noordermeer, and Dr. Emily E. Richards for useful feedback on the current literature. 

# References

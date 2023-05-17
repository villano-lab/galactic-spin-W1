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
  - name: Francis Vititoe 
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
&mdash; that is, these galaxies contain dark matter [@Rubin1978]. 

Participants in this workshop will have the opportunity to explore dark matter through scientific
literature-based [@Jimenez2003;@Karukes2015;@Richards2015;@Fraternali2011;@de_Naray2008] galactic rotation
curves both by using interactive programs and by editing Python code. This will give participants
an understanding of how physicists arrived at the idea of dark matter, showing them the difference
between curve fits with and without dark matter components. Understanding dark matter's
epistemological origins will help participants formulate their own opinions on the dark matter
debate.

## Materials

This project consists of several modules in the form of Jupyter notebooks [@Loizides2016]:

| **File Name**                                                    | Description                    |
| ---------------------------------------------------------------- | ------------------------------ |
| **01_DM_Rotation_Curve_Intro.ipynb**       | Animations and rotation curve plots demonstrating three types of rotational motion. |
|                                            |
| **02_Widget_NGC5533_DMonly.ipynb**         | Interactive widget to introduce dark matter. |
|                                            |
| **03_Measured_Data_Plotting.ipynb**        | Rotation curve plotting of measured velocities to visualize star and gas motions in a galaxy. |
|                                            |
| **04_Plotting_Rotation_Curves.ipynb**      | Plotting the rotation curves of galaxy components. |
|                                            |
| **05_Widget_NGC5533_All_Components.ipynb** | Interactive widget to visualize the components of the galaxy NGC 5533. |
|                                            |
| **06_Plotting_SPARC_Data.ipynb**           | Plotting the components of galactic rotation curves using the SPARC database of 175 galaxies. |
|                                            |
| **07_Bonus_Bulge_Rotation_Curve.ipynb**    | Constructing a rotation curve for the bulge component using empirically-derived parameters. |
|                                            |
| **08_Interactive_Fitting.ipynb**           | Interactive curve fitting.             |
|                                            |
| **09_Widget_SPARC_Galaxies.ipynb**         | Interactive widget to visualize the components of multiple galaxies using the SPARC database. |
|                                            |
| **10_Bonus_Black_Holes_as_DM.ipynb**       | Considering tiny black holes as dark matter candidates. |

# Statement of Need 

The primary goal of our project is to present rotation curve development and research in a versatile and approachable format 
for anyone to explore, learn from, and build upon. 
Rotation curves are a key empirical artifact through which dark matter can be observed and analyzed [@Rubin1978];
however, a thorough, start-to-finish description of the rotation curve building process is typically not given in scientific publications. Furthermore, software tools used in rotation curve literature are generally difficult for inexperienced users; 
for example, the GIPSY software package is very thorough but does not provide any introduction 
as it is intended for experienced users with a firm grasp on rotation curve components [@Gipsy1992]. 
Therefore, a rigorous yet accessible learning module is needed to provide an entry point 
for any individual interested in investigating the effect of dark matter in spiral galaxies. 
Our workshop is designed to present a convenient platform for developing basic rotation curves 
focused on introducing newcomers to the concepts necessary for understanding galactic rotation. 
This is achieved by leading users through hands-on computational activities, 
including building and plotting their own rotation curves. 

# Learning Objectives 

The learning objectives for these modules are:

1. Provide a working space where people can connect with current literature and identify as scientists.
2. Educate curious students or other individuals on the basic concepts of rotation curves, as related to the current problems and mysteries regarding dark matter in the universe.
3. Provide users with accessible activities relating to the basic principles of rotation curve composition. This includes:
     a. facilitating the introduction of rotation curve concepts via open-source code.
     b. interactive programs to provide users with practical and tangible approach of what producing rotation curves involves.
4. Understand data and models by interacting directly with equations and figures.

Most of the content provided in these modules has been presented and taught in previous
workshops/research symposiums (University of Colorado Denver: Data Science Symposium 2021 [@DataScienceSymposium2021],
Research and Creative Activities Symposium 2020 [@RaCAS2020], 2021 [@RaCAS2021], and 2022 [@RaCAS2022]) with feedback collected from
participants. We have chosen the activities for this module that proved most successful in terms
of education and sparking interest. 

# Delivery

The modules are designed to be presented to participants in numeric order as part of a workshop,
skipping those marked as "Bonus" as needed to fit the alloted time.  Participants are encouraged
to work together to complete the modules and compare their results to one another.  While working
through a module, the instructor(s) should be available to answer questions and check in on
participants' progress, but they should leave the bulk of the work to the participants themselves.
If any bonus modules are being skipped, the instructor(s) may wish to suggest them to participants
who find they are completing the content ahead of schedule.

All materials are designed to work on myBinder.org [@Binder], a website for hosting and interacting with jupyter notebooks. 
This is done to allow people to participate in the workshop without needing to install any software beforehand and 
is treated as the default delivery method. Participants who are experienced with python and 
already have a jupyter environmental installed may choose instead to run the modules locally.

# Story 

This project emerged from years of literary analysis and studying the reproducibility 
of rotation curve research. This journey impressed upon us a lack of clarity and accessibility for
newcomers in the world of rotation curves, not only in publications, but also in using software
and acquiring pre-existing data for rotation curve composition. The problems we encountered stemmed 
from the resources we found being very dense with technical language, 
focusing heavily on one or two components or even parameters, 
or assuming the reader has a certain level of familiarity with the subject prior to finding the resource in question. 
These traits are favorable for scientific journal content, but the lack of other types of content made it difficult to
find an entry point to the field. 
Our solution at the time was to dig into Noordermeer's paper on flattened Sérsic bulges [@Noordermeer2008], 
a paper that took us roughly a year to reproduce as we followed chains of references, corresponded with authors, and tried out rotation curve construction software in order to understand each rotation curve component. 
What we hope to accomplish is to provide others with necessary vocabulary and background knowledge before prompting them
to explore this kind of literature. Based on this experience and on feedback from our previous workshops and presentations
at research symposiums, 
we have developed our own software with a focus on improving accessibility and users' understanding
of the material by being clear, concise, and easily reproducible. 

# Acknowledgements

The authors would like to thank Dr. Martin Vogelaar at Kapteyn Astronomical Institute, Dr. Edo
Noordermeer, and Dr. Emily E. Richards for useful feedback on the current literature. 

# References

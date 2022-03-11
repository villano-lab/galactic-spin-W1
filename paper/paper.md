---
title: 'The Data Behind Dark Matter: Exploring Galactic Rotation' 
bibliography: references.bib
tags:
  - Python
  - Jupyter 
  - Binder
  - Workshop 
authors:
  - name: A.N. Villano
    affiliation: 1
    orcid: 0000-0002-3893-7259
  - name: Raphael Hatami
    affiliation: 1
  - name: Kitty C. Harris
    affiliation: 2
    orcid: 0000-0001-5406-8367
  - name: Judit Bergfalk 
    affiliation: 3
    orcid: 0000-0003-1662-0768 
  - name: F. Vititoe 
    affiliation: 4
affiliations:
 - name: Department of Physics, University of Colorado Denver, Denver CO 80217, USA
   index: 1
 - name: Integrated Sciences, University of Colorado Denver, Denver CO 80217, USA
   index: 2
 - name: Astrophysical & Planetary Sciences, University of Colorado Boulder, Boulder, CO 80309, USA 
   index: 3
 - name:  Department of Physics, University of Colorado Boulder, Boulder, CO 80309, USA 
   index: 4
date: 19 November 2021
nocite: '@*'
---

# Summary

By analyzing the rotational velocities of bodies in galaxies, physicists and astronomers have
found that there seems to be something missing in our understanding of these galaxies. One theory
is that there is some matter present in this galaxies which we cannot see because it doesn't
interact with light - that is, that these galaxies contain dark matter [@1978ApJ-225L-107R]. 

Participants in this workshop will have the opportunity to explore dark matter through scientific
literature-based [@Karukes_2015;@Richards_2015;@Fraternali_2011;@de_Naray_2008] galactic rotation
curves both by using interactive programs and by editing python code. This will give participants
an understanding of how physicists arrived at the idea of dark matter showing them the difference
between curve fits with and without galactic dark matter components. Understanding dark matter's
epistemological origins will help participants to formulate their own opinions on the dark matter
debate.

## Description of the software or learning module

In this interactive module, we include several programs in the form of Jupyter notebooks [@soton403913]:

| File Name | Short Description |
| --- | --- |
| 01_DM_Rotation_Curve_Intro.ipynb | Animations and rotation curve plots demonstrating three types of rotational motion. |
| 02_Widget_NGC5533_DMonly.ipynb | An interactive introduction to dark matter. | 
| 03_Measured_Data_Plotting.ipynb | Rotation curve plotting of measured velocities to visualize star and gas motions in a galaxy. | 
| 04_Plotting_Rotation_Curves.ipynb | Plotting the rotation curves of galaxy components. | 
| 05_Widget_NGC5533_All_Components.ipynb | Interactive widget to visualize the components of the galaxy NGC 5533. | 
| 06_Plotting_SPARC_Data.ipynb | Plotting the components of galactic rotation curves using the SPARC database of 175 galaxies. | 
| 07_Bonus_Bulge_Rotation_Curve.ipynb| Constructing a rotation curve for the bulge component using empirically-derived parameters. | 
| 08_Interactive_Fitting.ipynb | Interactive curve fitting. | 
| 09_Widget_SPARC_Galaxies.ipynb | Interactive widget to visualize the components of multiple galaxies using the SPARC database of 175 galaxies. | 
| 10_Bonus_Black_Holes_as_DM.ipynb | Considering tiny black holes as dark matter candidates. | 


# Statement of Need 

Rotation curves present one of the key empirical artifacts through which dark matter can be
observed and analyzed [@1978ApJ-225L-107R]. However, common understanding of rotation curves is lacking within
the field of astrophysics (source), especially in students and young academics in this field
(source). Furthermore, the software and tools used in actually developing rotation curves fall
into an even deeper level esotericism, exacerbated by various factors including poor accessibility
and outdated code (source: GIPSY) and lack of open-source resources available. We therefore have
created a Github repository consisting of Jupyter notebooks, data files, and notes comprising a
module of rotation curve learning. (OUR PROGRAM/MODULE NAME) is designed to present an
easily-accessible open-source platform for developing basic rotation curves in an annotated,
tutorial form for both academics and curious individuals with little or no experience with
rotation curves. With MODULE NAME, users can follow step-by-step instructions and explanations in
developing their own rotation curves. The primary goal of this project is to present
rotation-curve development and research in a versatile and approachable format for any individual
to explore, learn, and build upon.

# Learning Objectives 

	• For learning modules, describe the learning objectives, content, instructional design, and experience of use in teaching and learning situations.
The learning objectives for this educational module are:
	1. Educate curious students or other curious individuals on the basic concepts of rotation curves, as related to the current problems and mysteries regarding dark matter in the universe.
	2. Provide users with an accessible work and activities relating to the basic principles of rotation curve composition.
		a. This includes facilitating the introduction of rotation curve concepts via open-source, interactive programs to provide users with practical and tangible entrees of what rotation curve work really involves.
	3. Provide code-based concepts and examples of how users can become more involved in rotation curve work as they wish.

experience of use in teaching and learning situations.

(Raphael). Most of the content provided in this module has been presented and taught in previous workshops/research symposiums (Data Science Symposium, RaCAS Research Symposium 2020 and 2021), with feedback collected from participants. We have chosen for this module the activities which proved most successful in terms of education and sparking interest in participants. 

# Story 

(Raphael) The story of this project emerged out of several years of literary analysis, trial and error, and reproducibility studies our team conducted on rotation curve research. One principal feeling our team drew from this journey was a perceived lack of clarity and accessibility in the world of rotation curve research, including that pertaining to programs and data requisite for actually composing rotation curves. Examining the difficulties encountered in our team's research, we developed a simple resource (an educational module) that is clear, concise, with curve-composition activities that can be easily reproducible by the user. From our experience analyzing a number of rotation curve publications, in addition to correspondences with active researchers in the field, and our experience with currently-used curve-composing programs, we have begun to develop our own versions of such programs, making improvements which we believe will greatly improve accessibility and understanding. We also include content which we produced for previous workshops and research symposiums that further facilitate rotation curve learning. 

#Materials

* SPARC data import
* Interactive Measured Data Plotting: 
    * Plotting radial velocity measurements of multiple galaxies in a single plot to compare shapes of the curves
* Interactive Rotation Curve Plotting:
    * Choose a galaxy to plot
    * Build your own rotation curve using the galaxy's components. 
    * Calculate black hole component
    * Calculate dark matter component
    * See how dark matter component affects the total curve 
* Interactive Fitting:
    * Choose a galaxy to plot
    * Find scaling factors by fitting them to the measured data points
* Widget for Dark Matter parameters affecting the rotation curve in NGC 5533

Can make the same widget for the other galaxies too, if we want....

# Acknowledgements

The authors would like to thank...

# References

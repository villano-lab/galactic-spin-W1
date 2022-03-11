[![CC BY 4.0][cc-by-shield]][cc-by] 
[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/villano-lab/galactic-spin-W1/HEAD?labpath=binder)
[![Build Status](https://app.travis-ci.com/villano-lab/galactic-spin-W1.svg?branch=master)](https://app.travis-ci.com/github/villano-lab/galactic-spin-W1)
[![Documentation Status](https://readthedocs.org/projects/galactic-spin-W1/badge/?version=latest)](https://galactic-spin-W1.readthedocs.io/en/latest/?badge=latest)

# Galactic Spin W1: The Basics of Rotation Curves

The first workshop of our Rotation Curves series, which covers the basics of what rotation curves are,
how they are measured, and how we model them by breaking them into components.

![A poster, reading: CHOOSE: EMPIRICAL OR THEORETICAL? Astrophysicisists have developed models of the 'galactic bulge', the bright central component of spiral galaxies. Some of these models depend on empirical data sets as input, while more theoretical models only require empirically-derived parameters as input. Empirically-based models can be more accurate, however, purely-mathematical models can help us extrapolate, understand, and extend our knowledge beyond what can only be directly observed. When making your own bulge component, you'll get to choose your path: empirical or theoretical rotation curve models. CHOOSE: WHICH GALAXY? (End text) At the bottom of the poster, there are three labeled galaxy images. The first: NGC 5533, Constellation: BOOTES, DISTANCE: 177 MLY, SIZE: 50,000 LY. The second: NGC 891, CONSTELLATION: CAMELOPARDALIS, DISTANCE: 27.3 MLY, SIZE: 120,000 LY. The last: NGC 7814, CONSTELLATION: PEGASUS, DISTANCE: 40 MLY, SIZE: 80,000 LY.](binder/images/bulge.png)

You can find more detailed documentation about the code involved in this workshop [here](https://galactic-spin-W1.readthedocs.io/en/latest), 
including [a page guiding you through the workshop itself](https://galactic-spin-w1.readthedocs.io/en/feature-open-source-documents/01_Getting_Started.html).

## OUTLINE

| Module Name                                                                | Notebook Name                         | Description |
| -------------------------------------------------------------------------- | ------------------------------------- | ----------- |
| Understanding Rotation Curves                                              | 1_DM_Rotation_Curve_Intro.ipynb       | Rotation curves are a measure of how fast an object is moving at a certain distance from the center. Three kinds of rotation curves are discussed in this notebook: rigid-body rotation, planetary rotation, and galactic flat rotation. With the help of animations and rotation curve plots, the three cases are compared to gain a better understanding of galactic rotation curves. |
| Introduction to dark matter                                                | 2_Widget_NGC5533_DMonly.ipynb         | Mass in a galaxy can be visualized by plotting its rotation curve. A galactic rotation curve is a measure of how fast the stars and gas move in the galaxy at a certain distance from the center. According to Newton's law of gravitation, objects orbiting the center of gravity should depend on the mass enclosed in the system. However, the theoretical rotation curve of the measured visible matter does not agree with the measured velocities of matter in some of the spiral galaxies. Dark matter is introduced to account for the "missing matter". |
| Plotting measured velocities                                               | 3_Measured_Data_Plotting.ipynb        | Plotting measured velocities	3_Measured_Data_Plotting.ipynb	First step in understanding a rotation curve is to plot it. Comparing the rotation curves of multiple galaxies is a good exercise to visualize the motions of stars and gas. The shapes of curves reveal the mass distribution in a galaxy. For example, higher velocity measurements in the central region indicates a supermassive black hole at the center of that galaxy. A flat rotation curve suggests the presence of a dark matter halo. |
| Plotting the components of galactic rotation curves (4 galaxies)           | 4_Plotting_Rotation_Curves.ipynb      | Theoretical rotation curves are computed using the velocities of each component of the galaxy. The velocities of the bulge, disk and gas are calculated from luminosities, surface brightness profiles, surface density profiles or mass models but these calculations are beyond the scope of this workshop. For this reason, the velocities of the three components are imported into the notebook. On the contrary, the rotation curve of the central black hole (point-mass rotation curve - yet another type of rotation curve) and the dark matter halo can be easily produced. Adding all components, the total velocity of only light matter can then be compared to the total velocity that includes the dark matter component. Is it possible to fit both curves to the measured data? |
| Interactive widget to visualize the components of NGC 5533                  | 5_Widget_NGC5533_All_Components.ipynb | The rotation curves of each component in the galaxy NGC 5533 can be scaled up and down using the interactive widget in this notebook. The best combination of the scaling parameters results in a good fit to the measured data points. To characterize the goodness of the fit, the reduced chi-squared value is calculated. A value close to 1 is an indication of a good fit. |
| Plotting the components of galactic rotation curves (SPARC - 175 galaxies)  | 6_Plotting_SPARC_Data.ipynb           | The Spitzer Photometry & Accurate Rotation Curves (SPARC) database provides pre-calculated velocities of the bulge, disk and gas in 175 galaxies, as well as the measured velocity data points. Analogous to the 4_Plotting_Rotation_Curves.ipynb activity, this notebook also gives a guide to calculating the missing dark matter component and compares the rotation curve of luminous matter to the total velocity with the dark matter component. |
| Creating a rotation curve of the bulge component                            | 7_Bonus_Bulge_Rotation_Curve.ipynb    | The bulge is the most luminous, central component of a spiral galaxy that contains densely packed stars and gas. Although the rotation curve of the bulge can be derived from luminosity measurements, more theoretical models only utilize empirically-derived parameters such as the central surface brightness, the total luminosity of the bulge, the concentration parameter, and a characteristic radius. With the use of these parameters and calculus, the theoretical rotation of the bulge can be derived. |
| Interactive fitting of 4 galaxies                                           | 8_Interactive_Fitting.ipynb           | Fitting is a statistical method to scale parameters of a function until it closely resembles the curve of data points. Selecting a galaxy out of four options, these free parameters can be adjusted to vary. The results of the fitting is shown and explained in this notebook. |
| Interactive widget to visualize the components of multiple galaxies (SPARC) | 9_Widget_SPARC_Galaxies.ipynb         | After selecting a galaxy from the SPARC database of 175 galaxies, each component can be scaled up and down using the sliders of the interactive widget. Additionally, the fit parameters of the best fit to the measured velocities is revealed, along with an image of the selected galaxy from the NASA SkyView database. |
| Tiny black holes as dark matter candidates                                  | 10_Bonus_Black_Holes_as_DM.ipynb      | How many tiny black holes can account for the missing mass called dark matter? The interactive widgets and the visual representations of the number and mass of black holes give a good explanation to this question. For comparison, two spiral galaxies are investigated: NGC 5533 and NGC 7814. |

## CITATION

If you decide to use this code, or if you want to add a reference to it, please cite the latest archived version (soonTM),

## VERSION HISTORY

12.11.2021: Release of [version 0.0.1](https://github.com/villano-lab/galactic-spin-W1/releases/tag/v0.0.1)

## AUTHORS & CONTACT

The authors of *Galactic Spin W1* are A.N. Villano, Kitty C. Harris, Raphael Hatami, Judit Bergfalk, F. Vititoe, and Eric Lorimor.

For questions, support, bug reports, or other suggestions, please open an [issue](https://github.com/villano-lab/galactic-spin-W1/issues).

## LICENSE
This work is licensed under a
[Creative Commons Attribution 4.0 International License][cc-by].

[![CC BY 4.0][cc-by-image]][cc-by]

[cc-by]: http://creativecommons.org/licenses/by/4.0/
[cc-by-image]: https://i.creativecommons.org/l/by/4.0/88x31.png
[cc-by-shield]: https://img.shields.io/badge/License-CC%20BY%204.0-lightgrey.svg

Or see the [LICENSE file](https://github.com/villano-lab/galactic-spin-W1/blob/master/LICENSE).

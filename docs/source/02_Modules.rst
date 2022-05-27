==========
2. Modules
==========

--------
Overview
--------

Description of the learning modules


+------------------------------------------+------------------------------------------------------------+
| **File Name**                            | **Short Description**                                      | 
+==========================================+============================================================+
| `01_DM_Rotation_Curve_Intro.ipynb`       | Animations and rotation curve plots demonstrating          |
|                                          | three types of rotational motion.                          |
+------------------------------------------+------------------------------------------------------------+
| `02_Widget_NGC5533_DMonly.ipynb`         | An interactive introduction to Dark Matter.                | 
+------------------------------------------+------------------------------------------------------------+
| `03_Measured_Data_Plotting.ipynb`        | Rotation curve plotting of measured velocities             |
|                                          | to visualize star and gas motions in a galaxy.             |
+------------------------------------------+------------------------------------------------------------+
| `04_Plotting_Rotation_Curves.ipynb`      | Choose between 4 galaxies and plot the rotation            |
|                                          | curve of each component and their total velocity.          |
+------------------------------------------+------------------------------------------------------------+
| `05_Widget_NGC5533_All_Components.ipynb` | Interactive widget to adjust the components of             |
|                                          | the galaxy NGC 5533 and see how velocity changes.          |
+------------------------------------------+------------------------------------------------------------+
| `06_Plotting_SPARC_Data.ipynb`           | Plotting the components of galactic rotation curves        |
|                                          | using the SPARC database of 175 galaxies.                  |
+------------------------------------------+------------------------------------------------------------+
| `07_Bonus_Bulge_Rotation_Curve.ipynb`    | Calculate the theoretical rotation curve of the            |
|                                          | bulge component using empirically derived parameters.      |
+------------------------------------------+------------------------------------------------------------+
| `08_Interactive_Fitting.ipynb`           | Calculate the fitting parameters of the rotation           |
|                                          | curve to determine the amount of Dark Matter needed.       |
+------------------------------------------+------------------------------------------------------------+
| `09_Widget_SPARC_Galaxies.ipynb`         | Interactive widget to visualize the components of          |
|                                          | multiple galaxies using the SPARC database of 175 galaxies.|
+------------------------------------------+------------------------------------------------------------+
| `10_Bonus_Black_Holes_as_DM.ipynb`       | Considering tiny black holes as Dark Matter candidates.    |
+------------------------------------------+------------------------------------------------------------+

--------------------------
01_DM_Rotation_Curve_Intro
--------------------------

**Understanding rotation curves** 

.. image:: https://mybinder.org/badge_logo.svg
   :target: https://mybinder.org/v2/gh/villano-lab/galactic-spin-W1/master?labpath=binder/01_DM_Rotation_Curve_Intro.ipynb
`See module on GitHub <https://github.com/villano-lab/galactic-spin-W1/blob/master/binder/01_DM_Rotation_Curve_Intro.ipynb>`_

Rotation curves are a measure of how fast an object is moving at a certain distance from the center. Three kinds of rotation curves are discussed in this notebook: rigid-body rotation, planetary rotation, and galactic flat rotation. With the help of animations and rotation curve plots, the three cases are compared to gain a better understanding of galactic rotation curves.

------------------------
02_Widget_NGC5533_DMonly
------------------------

**Introduction to Dark Matter**

.. image:: https://mybinder.org/badge_logo.svg
   :target: https://mybinder.org/v2/gh/villano-lab/galactic-spin-W1/v1.0.0?labpath=binder%2F02_Widget_NGC5533_DMonly.ipynb
`See module on GitHub <https://github.com/villano-lab/galactic-spin-W1/blob/master/binder/02_Widget_NGC5533_DMonly.ipynb>`_

Mass in a galaxy can be visualized by plotting its rotation curve. A galactic rotation curve is a measure of how fast the stars and gas move in the galaxy at a certain distance from the center. According to Newton's law of gravitation, objects orbiting the center of gravity should depend on the mass enclosed in the system. However, the theoretical rotation curve of the measured visible matter does not agree with the measured velocities of matter in some of the spiral galaxies. Dark Matter is introduced to account for the "missing matter".

-------------------------
03_Measured_Data_Plotting
-------------------------

**Plotting measured velocities**

.. image:: https://mybinder.org/badge_logo.svg
   :target: https://mybinder.org/v2/gh/villano-lab/galactic-spin-W1/v1.0.0?labpath=binder%2F03_Measured_Data_Plotting.ipynb
`See module on GitHub <https://github.com/villano-lab/galactic-spin-W1/blob/master/binder/03_Measured_Data_Plotting.ipynb>`_

First step in understanding a rotation curve is to plot it. Comparing the rotation curves of multiple galaxies is a good exercise to visualize the motions of stars and gas. The shapes of curves reveal the mass distribution in a galaxy. For example, higher velocity measurements in the central region indicates a supermassive black hole at the center of that galaxy. A flat rotation curve suggests the presence of a Dark Matter halo. 

---------------------------
04_Plotting_Rotation_Curves
---------------------------

**Plotting the components of galactic rotation curves (4 galaxies)**

.. image:: https://mybinder.org/badge_logo.svg
   :target: https://mybinder.org/v2/gh/villano-lab/galactic-spin-W1/v1.0.0?labpath=binder%2F04_Plotting_Rotation_Curves.ipynb 
`See module on GitHub <https://github.com/villano-lab/galactic-spin-W1/blob/master/binder/04_Plotting_Rotation_Curves.ipynb>`_

Theoretical rotation curves are computed using the velocities of each component of the galaxy. The velocities of the bulge, disk and gas are calculated from luminosities, surface brightness profiles, surface density profiles or mass models but these calculations are beyond the scope of this workshop. For this reason, the velocities of the three components are imported into the notebook. On the contrary, the rotation curve of the central black hole (point-mass rotation curve - yet another type of rotation curve) and the dark matter halo can be easily produced. Adding all components, the total velocity of only light matter can then be compared to the total velocity that includes the dark matter component. Is it possible to fit both curves to the measured data?

--------------------------------
05_Widget_NGC5533_All_Components
--------------------------------

**Interactive widget to visualize the components of NGC 5533**

.. image:: https://mybinder.org/badge_logo.svg
   :target: https://mybinder.org/v2/gh/villano-lab/galactic-spin-W1/v1.0.0?labpath=binder%2F05_Widget_NGC5533_All_Components.ipynb 
`See module on GitHub <https://github.com/villano-lab/galactic-spin-W1/blob/master/binder/05_Widget_NGC5533_All_Components.ipynb>`_

The rotation curves of each component in the galaxy NGC 5533 can be scaled up and down using the interactive widget in this notebook. The best combination of the scaling parameters results in a good fit to the measured data points. To characterize the goodness of the fit, the reduced chi-squared value is calculated. A value close to 1 is an indication of a good fit.

----------------------
06_Plotting_SPARC_Data
----------------------

**Plotting the components of galactic rotation curves (SPARC - 175 galaxies)**

.. image:: https://mybinder.org/badge_logo.svg
   :target: https://mybinder.org/v2/gh/villano-lab/galactic-spin-W1/v1.0.0?labpath=binder%2F06_Plotting_SPARC_data.ipynb 
`See module on GitHub <https://github.com/villano-lab/galactic-spin-W1/blob/master/binder/06_Plotting_SPARC_data.ipynb>`_

The Spitzer Photometry & Accurate Rotation Curves (SPARC) database provides pre-calculated velocities of the bulge, disk and gas in 175 galaxies, as well as the measured velocity data points. Analogous to the 4_Plotting_Rotation_Curves.ipynb activity, this notebook also gives a guide to calculating the missing dark matter component and compares the rotation curve of luminous matter to the total velocity with the dark matter component. 

-----------------------------
07_Bonus_Bulge_Rotation_Curve
-----------------------------

**Creating a rotation curve of the bulge component**

.. image:: https://mybinder.org/badge_logo.svg
   :target: https://mybinder.org/v2/gh/villano-lab/galactic-spin-W1/v1.0.0?labpath=binder%2F07_Bonus_Bulge_Rotation_Curve.ipynb 
`See module on GitHub <https://github.com/villano-lab/galactic-spin-W1/blob/master/binder/07_Bonus_Bulge_Rotation_Curve.ipynb>`_

The bulge is the most luminous, central component of a spiral galaxy that contains densely packed stars and gas. Although the rotation curve of the bulge can be derived from luminosity measurements, more theoretical models only utilize empirically-derived parameters such as the central surface brightness, the total luminosity of the bulge, the concentration parameter, and a characteristic radius. With the use of these parameters and calculus, the theoretical rotation of the bulge can be derived. 

----------------------
08_Interactive_Fitting
----------------------

**Interactive fitting of 4 galaxies**

.. image:: https://mybinder.org/badge_logo.svg
   :target: https://mybinder.org/v2/gh/villano-lab/galactic-spin-W1/v1.0.0?labpath=binder%2F08_Interactive_Fitting.ipynb 
`See module on GitHub <https://github.com/villano-lab/galactic-spin-W1/blob/master/binder/08_Interactive_Fitting.ipynb>`_

Fitting is a statistical method to scale parameters of a function until it closely resembles the curve of data points. Selecting a galaxy out of four options, these free parameters can be adjusted to vary. The results of the fitting is shown and explained in this notebook. 

------------------------
09_Widget_SPARC_Galaxies
------------------------

**Interactive widget to visualize the components of multiple galaxies (SPARC)**

.. image:: https://mybinder.org/badge_logo.svg
   :target: https://mybinder.org/v2/gh/villano-lab/galactic-spin-W1/v1.0.0?labpath=binder%2F09_Widget_SPARC_Galaxies.ipynb 
`See module on GitHub <https://github.com/villano-lab/galactic-spin-W1/blob/master/binder/09_Widget_SPARC_Galaxies.ipynb>`_

After selecting a galaxy from the SPARC database of 175 galaxies, each component can be scaled up and down using the sliders of the interactive widget. Additionally, the fit parameters of the best fit to the measured velocities is revealed, along with an image of the selected galaxy from the NASA SkyView database. 

--------------------------
10_Bonus_Black_Holes_as_DM
--------------------------

**Tiny black holes as dark matter candidates**

.. image:: https://mybinder.org/badge_logo.svg
   :target: https://mybinder.org/v2/gh/villano-lab/galactic-spin-W1/v1.0.0?labpath=binder%2F10_Bonus_Black_Holes_as_DM.ipynb
`See module on GitHub <https://github.com/villano-lab/galactic-spin-W1/blob/master/binder/10_Bonus_Black_Holes_as_DM.ipynb>`_

How many tiny black holes can account for the missing mass called dark matter? The interactive widgets and the visual representations of the number and mass of black holes give a good explanation to this question. For comparison, two spiral galaxies are investigated: NGC 5533 and NGC 7814.
---
title: The Data Behind Dark Matter: Exploring Galactic Rotation 
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
affiliations:
 - name: Department of Physics, University of Colorado Denver, Denver CO 80217, USA
   index: 1
 - name: Integrated Sciences, University of Colorado Denver, Denver CO 80217, USA
   index: 2
 - name: Astrophysical & Planetary Sciences, University of Colorado Boulder, Boulder, CO 80309, USA 
   index: 2
date: 19 November 2021
nocite: '@*'
---

# Summary

By analyzing the rotational velocities of bodies in galaxies, physicists and astronomers have
found that there seems to be something missing in our understanding of these galaxies. One theory
is that there is some matter present in this galaxies which we cannot see because it doesn't
interact with light - that is, that these galaxies contain dark matter. Participants in this
workshop will have the opportunity to explore dark matter through galactic rotation curves both by
using interactive programs and by editing python code. This will give participants an
understanding of how physicists arrived at the idea of dark matter showing them the difference
between curve fits with and without dark matter. Understanding dark matter's epistemological
origins will help participants to formulate their own opinions on the dark matter debate.


# Models Used

When modeling deposits from neutron capture events, we want to look at the recoil of the nucleus
as a result of these cascades.  To determine how much energy is deposited, we must track how
much the nucleus slows down between steps of the cascade as well as *how* each state change
affects the nucleus' travel.  `nrCascadeSim` assumes a constant deceleration that results from the
nucleus colliding with other nearby nuclei.  This means that it must simulate, along with the
steps of the cascade, the time between each state &mdash; to calculate how much the nucleus slows
down &mdash; and the angle between the nucleus' momentum before a decay and the momentum boost
(gamma ray) resulting from the decay &mdash; to calculate the resulting momentum.  The time
between steps is simulated as an exponentially-decaying random variable based on the state's
half-life\footnote{It is most correct to use the half-life for the state given the state it will decay to. 
However, these are not generally well-known unless the branching ratios are well-known. 
If the ratios are well-known, then a correction can be made and incorporated into the input file.}, 
and the angle is simulated as having a uniform distribution on the surface of a sphere.
Cascade selection is weighted by isotope abundance and cross-section as well as the probability of
the energy level.  In+ existing levelfiles, energy levels are derived from 

The above process models the recoil energies, and the output gives both the total recoil energy
for a cascade as well as the energy per step.  For some applications, this may be the desired
output, or the user may already have a particular process they will use for converting this
energy to what they wish to measure.  However, we also include, for convenience, the ionization yield
and ionization energy of these recoils.  This ionization yield assumes the Lindhard
model:

$$
\begin{array}{rcl}
  Y & = & \frac{kg_{(\epsilon)}}{1+kg_{(\epsilon)}} \\
  g_{(\epsilon)} & = & a\epsilon^\gamma + b\epsilon^w + \epsilon \\
  \epsilon_{(E_r)} & = & 11.5E_r[keV]Z^{-7/3}
\end{array}
$$

Using the accepted value for Silicon ($k=0.143$) or Germanium ($k=0.159$), whichever is
appropriate; $a=3$; and $b=0.7$.

# Statement of Need

`nrCascadeSim` is a C/C++ package for generating a specified number of energy deposits resulting
from nuetron capture-induced nuclear recoils.  The energy levels and their lifetimes are
customizable, and multiple isotopes of the same element can be present within the simulation.
Pre-defined energy level files exist for silicon and germanium, which are constructed from the
data in .  Outputs include energy deposits at each step, total
kinetic energy deposits, and ionization energy deposits, making them useful for a variety of
applications, including nuclear recoil calibrations for dark matter direct detection or coherent
neutrino detection (CE$\mathrm{\nu}$NS).

# Example Use Case

Included in the repository is an example `test-example/Yields_and_Resolutions.ipynb` which users
can follow to ensure the code is running correctly.  This example both applies a yield model to
the individual energy deposits and applies variation intended to simulate the resolution of the
detector.  The yield and resolution models are described in more detail in the example notebook.
Figure \ref{LindvSor_fig} shows overlaid histograms of different combinations of analysis on the same
data file.

![An overlaid histogram showing an example use case in which points are generated and then multiple yield models and resolutions are applied.  In this example, the x-axis represents the ionization energy "yielded" by the cascade; this is effectively a way of noting what the detector reads out as opposed to what the pure kinetic energy of the cascade is.  The Lindhard yield  is output by \texttt{nrCascadeSim} as \texttt{Ei}; the Sorenson yield  is applied to the values from \texttt{delE}.Resolutions are applied by adding random values generated from a Gaussian distribution of fixed width to the energy yield.  The "Small Res (1/5)" histograms have Gaussians with 1/5 of the width of their counterparts.  The y-axis represents the normalized frequency of energy yields.\label{LindvSor_fig}](SorVsLin_fig.pdf)

# Acknowledgements

This material is based upon work supported by the U.S. Department of Energy, Office of Science, Office of High Energy Physics (HEP) under Award Number DE-SC0021364.

# References

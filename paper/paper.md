---
title: '`nrCascadeSim` - A simulation tool for nuclear recoil cascades resulting from neutron capture'
bibliography: references.bib
tags:
  - C++
  - Simulation
  - Nuclear Physics
authors:
  - name: A.N. Villano
    affiliation: 1
    orcid: 0000-0002-3893-7259
  - name: Kitty Harris
    affiliation: 1
    orcid: 0000-0001-5406-8367
  - name: Staci Brown
    affiliation: 2
affiliations:
 - name: Department of Physics, University of Colorado Denver, Denver CO 80217, USA
   index: 1
 - name: Department of Applied Mathematics & Statistics, University of New Mexico, Albuquerque NM 87131, USA
   index: 2
date: 12 December 2020
nocite: '@*'
---

# Summary

Neutron capture-induced nuclear recoils have emerged as an important
tool for detector calibrations in direct dark matter detection and coherent elastic neutrino-nucleus scattering (CE$\mathrm{\nu}$NS).

`nrCascadeSim` is a command-line tool for generating simulation data for energy deposits
resulting from neutron capture on pure materials. Presently, silicon, germanium, neon, and argon are
supported. While the software was developed for solid state detector calibration, it can be used
for any application which requires simulated neutron capture-induced nuclear recoil data.

A "cascade" occurs when a neutron becomes part of a nucleus.  The neutron can be captured to one
of many discrete energy levels, or states; if the energy level is nonzero (not the ground state),
then the state will eventually change so that it is zero.  This can happen either all at once or in
multiple steps &mdash; that is, the captured neutron may go from its state to the ground state, or
it may go to another state with lower energy that is not the ground state (provided that one
exists).  The cascade refers to the particular "path" of energy levels that a captured neutron
takes to get to the ground state from the neutron separation energy. Currently the code assumes
that the neutrons that enter the nuclear system have negligible (zero) kinetic energy; this is a
good approximation for thermal neutrons because 0.0254\ eV (the average kinetic energy of a
thermal neutron) is small compared to most nuclear recoil energy scales.

`nrCascadeSim` models many of these cascades at once and saves the energies along with other
useful data to a single file, the structure of which is outlined in Figure \ref{rootfile_fig}.

![An outline of the structure of a ROOT [@ROOT] output file named \texttt{file.root}. Everything is contained within a top-level key called \texttt{cascade}. Beneath \texttt{cascade} are several other keys, each pointing to an array. Each array element corresponds to one cascade; the same index will point to the same cascade across arrays. \texttt{n} notes the number of energy levels in the cascade. \texttt{cid} is short for "cascade ID" and refers to the row number of the levelfile which was used to generate the cascade, starting from zero. Each element of \texttt{Elev} is an array noting the energy levels used, given in eV. Similarly, \texttt{taus} notes the lifetimes of these states used, given in attoseconds. Both \texttt{Elev} and \texttt{taus} will have entries with a length of the corresponding value of n, so if \texttt{n[3]} is four then the lengths of \texttt{Elev[3]} and \texttt{taus[3]} will both be four. \texttt{delE} lists the energies deposited during the cascade in eV, and will always a length of one less than n. \texttt{I} calculates the ionization in terms of a number of charges, and \texttt{Ei} combines \texttt{I} with \texttt{delE} to list the ionization energy in eV. \texttt{time} describes the simulation-generated time that the neutron spent at each energy level, in attoseconds, and has a length corresponding to n. \texttt{Eg} provides gamma energies associated with each decay, in MeV, and has a length corresponding to one less than n. The gamma energies are not included in any of the other energy arrays. \label{rootfile_fig}](joss_fig.pdf)

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
the energy level.  In+ existing levelfiles, energy levels are derived from [@Ge] for germanium
and from [@Si] for Silicon.

The above process models the recoil energies, and the output gives both the total recoil energy
for a cascade as well as the energy per step.  For some applications, this may be the desired
output, or the user may already have a particular process they will use for converting this
energy to what they wish to measure.  However, we also include, for convenience, the ionization yield
and ionization energy of these recoils.  This ionization yield assumes the Lindhard
model[@lindhard]:

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
data in [@abundances] and [@nudat2].  Outputs include energy deposits at each step, total
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

![An overlaid histogram showing an example use case in which points are generated and then multiple yield models and resolutions are applied.  In this example, the x-axis represents the ionization energy "yielded" by the cascade; this is effectively a way of noting what the detector reads out as opposed to what the pure kinetic energy of the cascade is.  The Lindhard yield [@lindhard] is output by \texttt{nrCascadeSim} as \texttt{Ei}; the Sorenson yield [@sorensen] is applied to the values from \texttt{delE}.Resolutions are applied by adding random values generated from a Gaussian distribution of fixed width to the energy yield.  The "Small Res (1/5)" histograms have Gaussians with 1/5 of the width of their counterparts.  The y-axis represents the normalized frequency of energy yields.\label{LindvSor_fig}](SorVsLin_fig.pdf)

# Acknowledgements

This material is based upon work supported by the U.S. Department of Energy, Office of Science, Office of High Energy Physics (HEP) under Award Number DE-SC0021364.

# References

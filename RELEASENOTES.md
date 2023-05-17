## Release (v3.0.2) Date 17.05.23
Fix issues with documentation.
* Fixed a formatting issue causing version history to run together.
* Fixed conflicting license use -- now all code is under MIT license.
* Fixed error in a reference entry in `paper/references.bib`.

## Hotfix (v3.0.1) Date 07.04.23
Hotfix for restoring Binder functionality.
* The repository now successfully builds on Binder again.
* Fixed CircleCI builds only retrieving files from master for testing; they now retrieve files from the commit being tested.

## Release (v3.0.0) Date 28.07.22

Create documentation oriented toward developers and reorganize libraries.
* Documented library functions.
* Condensed small variables and library functions that are only called by other library functions into fewer, larger ones.
* Fixed readthedocs badges.

## Release (v2.0.1) Date 30.06.22

* Improved notebook 06 by restoring the ability to "Run all cells below" with a button press instead of JupyterLab UI.
* Added another CI to test an environment generated more similarly to how Binder generates its environments.

## Release (v2.0.0) Date 24.06.22

Condensed many library functions.
* PR #25 combined many independent but similar functions and libraries using components.py.
* PR #25 also fixed an issue in which the "Reveal!" button in notebook #5 sent its output to the log instead of printing normally.
* Fixed some inconsistencies in date format between different documents -- all dates should now be dd.mm.yy

## Release (v1.0.4) Date 17.06.22

Author name update.

## Release (v1.0.3) Date 17.06.22

Addressed author comments regarding current paper. 
* Citation changes
* Table format
* Clarifications/corrections

## Release (v1.0.2) Date 09.06.22

More thorough documentation, including updates to the paper.
* PR #19 added coverage data and badges.
* PR #20 created documentation for the notebooks to be used for the workshop.
* PR #21 completed a draft of the paper for author comments period.

## Release (v1.0.1) Date 27.04.22 

* Merged branches to fix missing documentation updates.

## Release (v1.0.0) Date 27.04.22 

DM Workshop 1 released for RaCAS 2022.
* Many files were renamed to improve the organization of the repository.
* Added Code of Conduct and Issue Templates.
* PR #3 fixed with the provided environment that sometimes caused sliders not to update anything.
* PR #4 added testing to Travis-CI
* PR #17 added documents such as CONTRIBUTING, PULL_REQUEST_TEMPLATE, and these RELEASENOTES.

## Release (v0.0.1) Date 08.01.22

* Initial release! Repository created.

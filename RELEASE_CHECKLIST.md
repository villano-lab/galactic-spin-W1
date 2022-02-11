#Release Checklist

Please be sure to do the following when making a release.

1. Update release page of readthedocs documentation `docs/source/07_Release_History.rst`.
2. Update citation page of readthedocs `docs/source/06_Citations.rst`.
3. Update the front-facing README in the Version history and citations sections.
5. Update release notes `RELEASENOTES.md`.
6. Update `docs/source/index.rst` if badges changed.
7. Restore fail-on-warning on .readthedocs.yaml if it was turned off.
8. Make a release on github.
9. Be sure the `stable` build of readthedocs points to the new release.
10. Be sure to create a version on readthedocs of the new release. 
11. If the version is a patch, deactivate the docs version for the previous patch of the same minor version. (Only one docs version for each minor version should be active at a time.)
12. Be sure codecov website is switched to default to master branch.
13. Update Zenodo.
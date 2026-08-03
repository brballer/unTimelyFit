# unTimelyFit
A search for Brown Dwarfs in the WISE satellite data

This search uses an independent proper motion (PM) fitting algorithm, unTimelyFit, than that used for the CatWISE catalog and an improved method for characterizing  proper motion likelihood, the survival function of the PM fit cumulative chi-squared distribution. It uses WISE IR data taken during the period 2010 through 2020. Candidates are initially assigned a quality classification using PM likelihood. The quality classification may be adjusted by the user after considering data from other instruments. Candidates are appended to a local catalog which allows for a later selection of objects for more detailed study.

The main Python module is doTile.py which launches processes to:
1) Identify the next "tile" (a patch of sky of ~1.6" x ~1.6") that has not been searched yet,
2) Download (public) unTimely and CatWISE source catalogs for this tile,
3) Download (public) DESI and Legacy Survey catalogs for this tile for later galaxy rejection,
4) Load (proprietary) unTimely calibration constants,
5) Re-process catWISE targets using processTile.py/unTimelyFit in multi-processor mode and merge results,
6) Write a list of star-like candidates having a high PM likelihood to a text file, and
7) Flag this tile as searched

The next, human driven, step is to individually confirm the star-like nature of each candidate in the text file by taking a poll of a dozen Vizier catalogs that use morphology for galaxy/star classification (dplan9.py). 

A candidate passing this requirement is visually inspected using a module that displays the positions, uncertainties and epochs of all unTimely objects within 3 pixels of the candidate position (plotUntimely.py). Objects are classified by a "recommended quality" variable (= 0, 1, 2, 3) that is an assessment of motion using the differential survival of the cumulative probability density function. A quality of 3 indicates a high differential probability that a moving hypothesis is correct. Objects having a low differential probability are most likely not moving despite some having a high PM significance for the moving hypothesis.  After reviewing data from other sources, e.g. the Legacy Survey Sky Browser, the user is prompted to accept or modify the recommended quality. The unTimelyFit results and the quality classification are stored in a local catalog if the quality is deemed significant.

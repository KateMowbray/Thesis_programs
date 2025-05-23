Changes to James Threlfall's particle orbit code allowing for work with collapsing magnetic trap (CMT) models

This repository contains each of the files that I have adjusted or written when working with the particle orbit code found at https://github.com/jwt104/party_orb/tree/master

The most major changes are the splitting of the CMT fields module into a 2D & 2.5D CMT module and a 3D CMT module along with adjustments made to make the new module CMTfields2half_mod compatible with the remaining code.

Other changes include the measurements of particular energisation terms in r_rkdrive_mod.f90 and the addition of missing terms in the db/dt term in r_derivs_mod.f90.

The remaining programs concern data visualisation, with scattercol4.pro giving the scatter plots showing the ratio of final to initial energies of particles against the ratio of final to initial field strength and the programs in the IDL folder being used to sketch 2D and 3D field lines for CMT models.

Visualisation for the breaking jet model is done with manyplots.pro and abenergies.pro, which show the energy profiles of a large number of particles and abscatter, which is similar to scattercol4, but grays out particles that escape early.

These files should be used with the code in the above link (which also contains a manual for the code), with files from this repository being used when they have the same name as files in that repository. The README for the linked repository had been included here as jtREADME.md.

Additional files describe the estimator used in the final full chapter in 2D (fermiguessfor.f90) and 3D (fermiguess3D.f90). The programs estaccuracytest.pro and estaccuracyplots.pro are used to quantify the errors for a large data set and to plot the locations associated with these errors, while fermiguessforcharts is used to sompare the energies of smaller data sets graphically.

The README gives further details on how to obtain the results presented in the thesis.

# Generate crystal deformation with Langragian strain for VASP calculations
To generate input files for VASP. It uses the lagrangian strain to strain the system with user specified deformation constant. This value should be in the harmonic regime otherwise this approximation will fail.

NB: This code is the modification of the script provided with theexciting code. All rights belongs to the original author,
if there is an error in the code or bug please contact the exciting team.

**This script run with python 3 or higher version.**
[![Maintenance](https://img.shields.io/badge/Maintained%3F-yes-green.svg)](https://GitHub.com/Naereen/StrapDown.js/graphs/commit-activity)
![pypi](https://img.shields.io/pypi/v/pybadges.svg)
![versions](https://img.shields.io/pypi/pyversions/pybadges.svg)
[![GPLv3 license](https://img.shields.io/badge/License-GPLv3-blue.svg)](http://perso.crans.org/besson/LICENSE.html)

I have modified the script to read CONTCAR optimized file (already minimized with IBRION=2, ISIF=3) and then print out the deformed POSCAR files with various strain.
Just follow the recipe on the screen and it will guide you to the rest of the process. For more information please visit exciting website (Elastic code).
The literature behind can be found at the exciting website.


ref: [exciting](http://exciting-code.org/nitrogen-energy-vs-strain-calculations)
ref: [materialsproject](https://wiki.materialsproject.org/Elasticity_calculations)

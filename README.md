# Crystal structure deformation using Langragian strain for VASP input files
**This script compiles with python 3 or higher version.**
[![Maintenance](https://img.shields.io/badge/Maintained%3F-yes-green.svg)](https://GitHub.com/Naereen/StrapDown.js/graphs/commit-activity)
![pypi](https://img.shields.io/pypi/v/pybadges.svg)
![versions](https://img.shields.io/pypi/pyversions/pybadges.svg)
[![GPLv3 license](https://img.shields.io/badge/License-GPLv3-blue.svg)](http://perso.crans.org/besson/LICENSE.html)

ENERGY Vs STRAIN APPROACH: Two scripts are provided for different deformation scheme. Please read the documentation or 
the details on the screen.

To generate input files for VASP. It uses the lagrangian strain to strain the system with user specified deformation constant. This value should be in the harmonic regime otherwise this approximation will fail. Large deformation leads to phase transitions.

NB: This code is the modification of the script provided with the exciting code. All rights belongs to the original author,
if there is an error in the code or bug please contact the exciting team.

```
**Deformation 1:
I have modified the script to read CONTCAR optimized file (already minimized with IBRION=2, ISIF=3) and then print out the deformed POSCAR files with various strain.
Just follow the recipe on the screen and it will guide you to the rest of the process. For more information please visit exciting website (Elastic code).
The literature behind can be found at the exciting website.
```

```
**Deformation 2:
Cubic_strain.py script uses another deformation matrix technique. Please note it only pertians to 
cubic system. It can be extended to include other crystal deformation matrix.
```

[1] ref: [exciting](http://exciting-code.org/nitrogen-energy-vs-strain-calculations)

[2] ref: [materialsproject](https://wiki.materialsproject.org/Elasticity_calculations)

**Another interesting is: STRESS Vs STRAIN Appraoch  

[3] ref: [ELASTIC](https://elastic.readthedocs.io/en/stable/index.html) 

[L.D. Landau, E.M. Lifszyc, Theory of elasticity, Elsevier (1986) ; ISBN: 075062633X, 9780750626330]

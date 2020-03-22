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
I have modified the script to read CONTCAR optimized file (already minimized with IBRION=2, ISIF=3) 
and then print out the deformed POSCAR files with various strain.
Just follow the recipe on the screen and it will guide you to the rest of the process. 
For more information please visit exciting website (Elastic code). The literature behind can be found at the exciting website.
	print ("                        |η[0]    η[5]/2  η[4]/2| ")
	print ("                    η = |η[5]/2  η[1]    η[3]/2| ")
	print ("                        |η[4]/2  η[3]/2  η[2]  | ")  
	print ("                               D' = I + η") 
```

```
**Deformation 2:
Cubic_strain.py script uses another deformation matrix technique. Please note it only pertians to 
cubic system. It can be extended to include other crystal deformation matrix.
```

[1] ref: [exciting](http://exciting-code.org/nitrogen-energy-vs-strain-calculations)

[2] ref: [materialsproject](https://wiki.materialsproject.org/Elasticity_calculations)

**Another interesting appraoch is: STRESS Vs STRAIN **

[3] ref: [ELASTIC](https://elastic.readthedocs.io/en/stable/index.html) 

[L.D. Landau, E.M. Lifszyc, Theory of elasticity, Elsevier (1986) ; ISBN: 075062633X, 9780750626330]

𝗢𝗻𝗲 𝗳𝗶𝗲𝗹𝗱 𝗼𝗳 𝘄𝗼𝗿𝗸 𝗶𝗻 𝘄𝗵𝗶𝗰𝗵 𝘁𝗵𝗲𝗿𝗲 𝗵𝗮𝘀 𝗯𝗲𝗲𝗻 𝘁𝗼𝗼 𝗺𝘂𝗰𝗵 𝘀𝗽𝗲𝗰𝘂𝗹𝗮𝘁𝗶𝗼𝗻 𝗶𝘀 𝗰𝗼𝘀𝗺𝗼𝗹𝗼𝗴𝘆. 𝗧𝗵𝗲𝗿𝗲 𝗮𝗿𝗲 𝘃𝗲𝗿𝘆 𝗳𝗲𝘄 𝗵𝗮𝗿𝗱 𝗳𝗮𝗰𝘁𝘀 𝘁𝗼 𝗴𝗼 𝗼𝗻, 𝗯𝘂𝘁 𝘁𝗵𝗲𝗼𝗿𝗲𝘁𝗶𝗰𝗮𝗹 𝘄𝗼𝗿𝗸𝗲𝗿𝘀 𝗵𝗮𝘃𝗲 𝗯𝗲𝗲𝗻 𝗯𝘂𝘀𝘆 𝗰𝗼𝗻𝘀𝘁𝗿𝘂𝗰𝘁𝗶𝗻𝗴 𝘃𝗮𝗿𝗶𝗼𝘂𝘀 𝗺𝗼𝗱𝗲𝗹𝘀 𝗳𝗼𝗿 𝘁𝗵𝗲 𝗨𝗻𝗶𝘃𝗲𝗿𝘀𝗲, 𝗯𝗮𝘀𝗲𝗱 𝗼𝗻 𝗮𝗻𝘆 𝗮𝘀𝘀𝘂𝗺𝗽𝘁𝗶𝗼𝗻𝘀 𝘁𝗵𝗮𝘁 𝘁𝗵𝗲𝘆 𝗳𝗮𝗻𝗰𝘆. 𝗧𝗵𝗲𝘀𝗲 𝗺𝗼𝗱𝗲𝗹𝘀 𝗮𝗿𝗲 𝗽𝗿𝗼𝗯𝗮𝗯𝗹𝘆 𝗮𝗹𝗹 𝘄𝗿𝗼𝗻𝗴. 𝗜𝘁 𝗶𝘀 𝘂𝘀𝘂𝗮𝗹𝗹𝘆 𝗮𝘀𝘀𝘂𝗺𝗲𝗱 𝘁𝗵𝗮𝘁 𝘁𝗵𝗲 𝗹𝗮𝘄𝘀 𝗼𝗳 𝗻𝗮𝘁𝘂𝗿𝗲 𝗵𝗮𝘃𝗲 𝗮𝗹𝘄𝗮𝘆𝘀 𝗯𝗲𝗲𝗻 𝘁𝗵𝗲 𝘀𝗮𝗺𝗲 𝗮𝘀 𝘁𝗵𝗲𝘆 𝗮𝗿𝗲 𝗻𝗼𝘄. 𝗧𝗵𝗲𝗿𝗲 𝗶𝘀 𝗻𝗼 𝗷𝘂𝘀𝘁𝗶𝗳𝗶𝗰𝗮𝘁𝗶𝗼𝗻 𝗳𝗼𝗿 𝘁𝗵𝗶𝘀. 𝗧𝗵𝗲 𝗹𝗮𝘄𝘀 𝗺𝗮𝘆 𝗯𝗲 𝗰𝗵𝗮𝗻𝗴𝗶𝗻𝗴, 𝗮𝗻𝗱 𝗶𝗻 𝗽𝗮𝗿𝘁𝗶𝗰𝘂𝗹𝗮𝗿, 𝗾𝘂𝗮𝗻𝘁𝗶𝘁𝗶𝗲𝘀 𝘁𝗵𝗮𝘁 𝗮𝗿𝗲 𝗰𝗼𝗻𝘀𝗶𝗱𝗲𝗿𝗲𝗱 𝘁𝗼 𝗯𝗲 𝗰𝗼𝗻𝘀𝘁𝗮𝗻𝘁𝘀 𝗼𝗳 𝗻𝗮𝘁𝘂𝗿𝗲 𝗺𝗮𝘆 𝗯𝗲 𝘃𝗮𝗿𝘆𝗶𝗻𝗴 𝘄𝗶𝘁𝗵 𝗰𝗼𝘀𝗺𝗼𝗹𝗼𝗴𝗶𝗰𝗮𝗹 𝘁𝗶𝗺𝗲. 𝗦𝘂𝗰𝗵 𝘃𝗮𝗿𝗶𝗮𝘁𝗶𝗼𝗻𝘀 𝘄𝗼𝘂𝗹𝗱 𝗰𝗼𝗺𝗽𝗹𝗲𝘁𝗲𝗹𝘆 𝘂𝗽𝘀𝗲𝘁 𝘁𝗵𝗲 𝗺𝗼𝗱𝗲𝗹 𝗺𝗮𝗸𝗲𝗿𝘀." 𝗗𝗶𝗿𝗮𝗰, 𝗣𝗮𝘂𝗹. 𝗢𝗻 𝗺𝗲𝘁𝗵𝗼𝗱𝘀 𝗶𝗻 𝘁𝗵𝗲𝗼𝗿𝗲𝘁𝗶𝗰𝗮𝗹 𝗽𝗵𝘆𝘀𝗶𝗰𝘀. (𝗧𝗿𝗶𝗲𝘀𝘁𝗲. 𝗝𝘂𝗻𝗲 𝟭𝟵𝟲𝟴 .) 

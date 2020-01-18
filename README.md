# Generate Langragian strain crystals for VASP calculation
Generating input files for VASP. Here it uses the lagrangian strain to strain the system with users specified deformation constant. This value should be in the harmonic
approximation otherwise it will fail.

NB: This code is the modification of the script provided with theexciting code. All rights belongs to the original author,
if there is an error in the code or bug please contact the exciting team.

**This script run with python 3 or higher version.**
I have modified the script to read CONTCAR optimized file (already minimized with IBRION=2, ISIF=3) and then print out the deformed POSCAR files with various strain.
Just follow the recipe on the screen and it will guide you to the rest of the process. For more information please visit exciting website (Elastic code).
The literature behind can be found at the exciting website.

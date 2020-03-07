#!/usr/bin/env python3

'''
A simple script to extract stress tensor from the vasp xml file and 
compute the pressure and correction to the stress tensor.

This script uses ase module to exttract data instead of writing the code
from scratch.
https://wiki.fysik.dtu.dk/ase/_modules/ase/io/vasp.html#read_vasp_xml
https://wiki.fysik.dtu.dk/ase/_modules/ase/atoms.html#Atoms.get_forces

'''
#-------------------------------------------------------------------------------#
import numpy as np
import ase.units
#from ase import Atoms
import ase.io
import ase
from numpy import dot, diag, ones, reshape, linspace, array, mean
from math import acos, pi, cos, sin, sqrt
#-------------------------------------------------------------------------------#
I = np.identity(3)
file = ase.io.read(filename='vasprun.xml',index=-1)
print ("From vasprun.xml file stress tensor is given in matrix notation:")
print ("Volume -> {:6.4f}".format(file.get_volume() ))

#print ("Force ->{}".format(file.get_forces() ))
print ("Lattice vectors::\n", end="")
print ("{}".format(np.mat(file.get_cell()) ), end="\n")

stress = file.get_stress(voigt=False) # this can be switched
print ("Stress tensor in Matrix notation:\n",stress)
p = -(stress.trace() )/3 # Hydrostatic Pressure
hyd_stress = (stress.trace() )/3 # Hydrostatic Stress

print ("pressure ->", p)
print ("Stress - pressure:\n",stress - dot(p,I) )
print ("="*80)
#-------------------------------------------------------------------------------#
'''
# This is the manual demonstration of how ase reads vasprun.xml file
# One thing to remember in vasp stress units are in kb which can be 
# converted to GPa by 0.1*GPa

  param :: pressure = - stress
'''
#-------------------------------------------------------------------------------#
stress1 = np.zeros((3, 3))
sblocks = [
   [ -0.13697103,       0.05274138,       0.24341161],
   [  0.05274208,      -0.65966872,       0.33726817],
   [  0.24341132,       0.33726859,      -0.62738792] ]
s = np.array(sblocks)	 
if s is not None:
	for i, v in enumerate(s):
		stress1[i] = np.array([float(val) for val in v[:]])
	stress1 *= -0.1 * ase.units.GPa	
	stress1 = stress1.reshape(9)[[0, 4, 8, 5, 2, 1]]
	print(stress1)
	p = - mean(stress1[:3]) # Hydrostatic Pressure	
	print ("pressure =>", p)
	print ("Stress - pressure =>\n", stress1 - p )
	

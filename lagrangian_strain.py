#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#_______________________________________________________________________________

from   numpy import *
import subprocess, shutil
import os.path
import numpy as np
import math, sys, os
import matplotlib.pyplot as plt
from colorama import Fore, Back, Style, init
init(autoreset=True)
#_______________________________________________________________________________
CRED = '\033[91m';CEND = '\033[0m'

def Elastic_strain():
	print ('\n{:_^80s}'.format("Documentation"))
	print ("__| {:10s}:: {:40s}".format("Author","Asif Iqbal"))
	print ("__| {:10s}:: {:40s}".format("DATED","06/03/2020"))
	print ("__| {:10s}:: {:40s}".format("USAGE","python3 sys.argv[0] <directory_name>"))	
	print (CRED+"__| [Please visit] http://exciting-code.org/nitrogen-energy-vs-strain-calculations" + CEND)
	print (CRED+"__| \u03B7  = \u03B5 + 1/2*\u03B5**2 :: Deformation scheme" + CEND)
	print (CRED+"__| r' = (I + \u03B5) * r :: Original to new coordinates axis" + CEND)
	print (CRED+"__| This script doesn't generae deformation based on Crystal symmetry rather" + CEND)
	print (CRED+"__| you can specify which deformation you want to investigate." + CEND)
	print (CRED+"__| Deformation scheme uses Voigt notation." + CEND)
	print ("                        |η[0]    η[5]/2  η[4]/2| ")
	print ("                    η = |η[5]/2  η[1]    η[3]/2| ")
	print ("                        |η[4]/2  η[3]/2  η[2]  | ")  
	print ("                               D' = I + η") 
	print ('{:_^80s}\n'.format("END of Documentation"))
	
	if (str(os.path.exists('CONTCAR'))=='False'): 
		sys.exit("ERROR: Input file CONTCAR not found!\n")
	
	maximum_strain = float( input("Enter maximum Lagrangian strain [1-10%] >>>> ") )
	if (1 < maximum_strain or maximum_strain < 0): 
		sys.exit("ERROR: Maximum Lagrangian strain is out of range [0-1]!\n")
	strain_points = int( input("Enter # of strain values (odd preferably) >>>> ") )
	tmp = int ( floor(strain_points/2) )
	print("The deformation range is [-{},{}]".format(tmp, tmp) )
	if (3 > strain_points or strain_points > 99): 
		sys.exit("ERROR: Number of strain values is out of range [3-99]!\n")
	
	print(Style.RESET_ALL)	
	print (Back.GREEN  + "------------------------------------------------------------------------" )
	print (Back.GREEN  + " List of deformation codes for strains in Voigt notation"                 )
	print (Back.GREEN  + "------------------------------------------------------------------------" )
	print (Back.YELLOW + "  0 =>  ( η, η, η, 0, 0, 0)  | volume strain "           )
	print (Back.YELLOW + "  1 =>  ( η, 0, 0, 0, 0, 0)  | linear strain along x "   )
	print (Back.GREEN  + "  2 =>  ( 0, η, 0, 0, 0, 0)  | linear strain along y "   )
	print (Back.GREEN  + "  3 =>  ( 0, 0, η, 0, 0, 0)  | linear strain along z "   )
	print (Back.GREEN  + "  4 =>  ( 0, 0, 0, η, 0, 0)  | yz shear strain"          )
	print (Back.GREEN  + "  5 =>  ( 0, 0, 0, 0, η, 0)  | xz shear strain"          )
	print (Back.GREEN  + "  6 =>  ( 0, 0, 0, 0, 0, η)  | xy shear strain"          )
	print (Back.YELLOW + "  7 =>  ( 0, 0, 0, η, η, η)  | shear strain along (111)" )
	print (Back.GREEN  + "  8 =>  ( η, η, 0, 0, 0, 0)  | xy in-plane strain "      )
	print (Back.GREEN  + "  9 =>  ( η,-η, 0, 0, 0, 0)  | xy in-plane shear strain" )
	print (Back.GREEN  + " 10 =>  ( η, η, η, η, η, η)  | global strain"            )
	print (Back.GREEN  + " 11 =>  ( η, 0, 0, η, 0, 0)  | mixed strain"             )
	print (Back.GREEN  + " 12 =>  ( η, 0, 0, 0, η, 0)  | mixed strain"             )
	print (Back.GREEN  + " 13 =>  ( η, 0, 0, 0, 0, η)  | mixed strain"             )
	print (Back.GREEN  + " 14 =>  ( η, η, 0, η, 0, 0)  | mixed strain"             )
	print (Back.GREEN  + "------------------------------------------------------------------------" )
	print(Style.RESET_ALL)	
	#-------------------------------------------------------------------------------		
	deformation_code = int(input("\nEnter deformation code >>>> "))
	if (0 > deformation_code or deformation_code > 14): 
		sys.exit("ERROR: Deformation code is out of range [0-14]!\n")
	
	if (deformation_code == 0 ): dc='EEE000'
	if (deformation_code == 1 ): dc='E00000'
	if (deformation_code == 2 ): dc='0E0000'
	if (deformation_code == 3 ): dc='00E000'
	if (deformation_code == 4 ): dc='000E00'
	if (deformation_code == 5 ): dc='0000E0'
	if (deformation_code == 6 ): dc='00000E'
	if (deformation_code == 7 ): dc='000EEE'
	if (deformation_code == 8 ): dc='EE0000'
	if (deformation_code == 9 ): dc='Ee0000'
	if (deformation_code == 10): dc='EEEEEE'
	if (deformation_code == 11): dc='E00E00'
	if (deformation_code == 12): dc='E000E0'
	if (deformation_code == 13): dc='E0000E'
	if (deformation_code == 14): dc='EE0E00'
	
	#-------------------------------------------------------------------------------
	file1 = open("CONTCAR",'r')
	line1 = file1.readlines()		
	file1.close()
	for i in line1:
		if ("Direct" or "direct" or "d" or "D") in i:
			PP=line1.index(i)
	#-------------------------------------------------------------------------------
	
	input_obj 	= open("CONTCAR","r")
	
	firstline   = input_obj.readline() # IGNORE first line comment
	scale 	= float(input_obj.readline())  # scale
	Latvec1 	= input_obj.readline()
	Latvec2 	= input_obj.readline()            
	Latvec3 	= input_obj.readline()            
	elementtype	= input_obj.readline().split()
	if (str.isdigit(elementtype[0])):
		sys.exit("VASP 4.X POSCAR detected. Please add the atom types")
	atom_number = input_obj.readline()
	Coordtype	= input_obj.readline()
	nat 		= atom_number.split()
	nat 		= [int(i) for i in nat]
	print ("Number of atoms in the cell:: {} ".format( sum(nat)) )
	
	input_obj.close()
	#-------------------------------------------------------------------------------
	a=[]; b=[]; c=[];
	Latvec1		= Latvec1.split()
	Latvec2		= Latvec2.split()
	Latvec3		= Latvec3.split()	
	for ai in Latvec1: 	a.append(float(ai))
	for bi in Latvec2: 	b.append(float(bi))
	for ci in Latvec3: 	c.append(float(ci))
	
	xml_basevect = np.array([a] + [b] + [c])		
	#print ("{}".format(xml_basevect),end="\n" )
	axis_matrix = np.array(xml_basevect) 
	determinant = np.linalg.det(axis_matrix)
	volume = np.abs(determinant*scale**3)
	print("Equilibrium volume of cell:: {} ".format(volume) )
	#-------------------------------------------------------------------------------
	work_directory = 'workdir'
	if (len(sys.argv) > 1): work_directory = sys.argv[1]
	if (os.path.exists(work_directory)): shutil.rmtree(work_directory)
	os.mkdir(work_directory)
	os.chdir(work_directory)
	
	output_info = open('INFO-elastic-constants',"w")
	
	output_info.write("\nMaximum Lagrangian strain       = {}".format( maximum_strain ))
	output_info.write("\nNumber of strain values         = {}".format(strain_points))
	output_info.write("\nVolume of equilibrium unit cell = {} [A]^3".format(volume))
	output_info.write("\nDeformation code                = {}".format(deformation_code))
	output_info.write("\nDeformation label               = {}".format(dc, "\n"))
	
	output_info.close()
	
	#-------------------------------------------------------------------------------
	
	delta=strain_points-1 ;# print (delta)
	convert=1
	eta_step=2*maximum_strain/delta
	#print(eta_step)
	#-------------------------------------------------------------------------------
	t = 1; tmp=-tmp;
	print ("{:12s} {:12.8s} {:14.8s} {:14.8s}".format("", "Vol_cell", "Vol_D'", "V/V_D'" ))
	for i in range(0,strain_points):
		eta=i*eta_step-maximum_strain*convert
		#print (eta)
		if (i+1 < 10): strainfile = 'strain-'+str(i+1).zfill(2)
		else: strainfile = 'strain-'+str(i+1).zfill(2)
		output_str = open(strainfile,"w")
		output_str.write( "{:10.5f}\n".format(eta) )
		output_str.close()
	
		if (abs(eta) < 0): eta=0
		ep=eta 
		if (eta < 0.0): em=abs(eta)
		else: em=-eta
		
	#-------------------------------------------------------------------------------
	
		e=[]
		for j in range(6):
			ev=0
			if  (dc[j:j+1] == 'E' ): ev=ep; #print (ev)
			elif(dc[j:j+1] == 'e' ): ev=em
			elif(dc[j:j+1] == '0' ): ev=0
			else: print ("==> "), dc; sys.exit("ERROR: deformation code not allowed!") 
			e.append(ev) 
		e = np.array(e)
		#print(numpy.array(e))	
	#-------------------------------------------------------------------------------
	
		eta_matrix=np.mat( [
		[ e[0]  , e[5]/2, e[4]/2], 
		[ e[5]/2, e[1]  , e[3]/2], 
		[ e[4]/2, e[3]/2, e[2]  ] ], dtype=float32 )
		one_matrix=np.identity(3)
	
	#-------------------------------------------------------------------------------
					
		norma=1 ; inorma=0 ; eps_matrix=eta_matrix
		if (np.linalg.norm(eta_matrix) > 0.7):sys.exit("ERROR: too large deformation!") 
	
		while ( norma > 1.e-10 ):
			x=eta_matrix - (1/2) * np.dot(eps_matrix,eps_matrix)
			norma=np.linalg.norm(x-eps_matrix)     
			#print(norma)		
			eps_matrix=x
			inorma=inorma+1
	
		def_matrix=one_matrix+eps_matrix	
		new_axis_matrix=np.transpose(np.dot(def_matrix,np.transpose(axis_matrix)))
		nam=np.mat( new_axis_matrix, dtype=float32  )
		V = np.linalg.det(new_axis_matrix)
		V_def = np.linalg.det(def_matrix)
		print ("{:02d}({:2d}) => {:10.6f} {:10.6f} {:14.6f}".format(t, tmp, abs(V), abs(V_def), V/V_def ) ) 

	#-------------------------------------------------------------------------------
	
		if (i+1 < 10): 
			outputfile = 'POSCAR-'+str(i+1).zfill(2)
		else: 
			outputfile = 'POSCAR-'+str(i+1).zfill(2)
		output_obj = open(outputfile,"w")
		output_obj.write(firstline)
		output_obj.write("{:10.8f}\n".format((scale)))
	
		for j in range(3):
			output_obj.write("{:22.16f} {:22.16f} {:22.16f}\n".format( (nam[j,0]), (nam[j,1]), (nam[j,2]) )  ) 
			
		for j in elementtype:
			output_obj.write("\t" +  j)
		output_obj.write("\n" )
		
		for j in atom_number:
			output_obj.write( "{}".format(j) )	
	
		for i in range(len(line1)-PP):
			output_obj.write(line1[PP+i] )	
		output_obj.close()
		
		t+=1; tmp +=1
#-------------------------------------------------------------------------------
	os.chdir('../')
	print ("\n")

if __name__ == "__main__":	
	Elastic_strain()
	
#plt.imshow(numpy.log(numpy.abs(numpy.fft.fftn(nam))**2))
#plt.show()

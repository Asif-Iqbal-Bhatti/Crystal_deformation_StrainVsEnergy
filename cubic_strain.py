#!/usr/bin/env python3
#_______________________________________________________________________________
from   lxml  import etree
from   sys   import stdin
from   numpy import *
import subprocess, os, sys
import os.path, shutil
import numpy as np 
import math, ase
#_______________________________________________________________________________
CRED = '\033[91m';CEND = '\033[0m'
CYEL = '\033[33m'; CEND = '\033[0m'
CPIN = '\033[44m';
#os.system('ase -T convert -i vasp CONTCAR -o exciting -f input.xml')
if (str(os.path.exists('input.xml'))=='False'): 
    sys.exit("ERROR: Input file input.xml not found!\n")

print ( "="*80 )
print ("__| {:10s}:: {:40s}".format("Author","Asif Iqbal"))
print ("__| {:10s}:: {:40s}".format("DATED","05/03/2020"))
print ("__| {:10s}:: {:40s}".format("USAGE","python3 sys.argv[0] <directory_name>"))
print ('{:_^80s}'.format("Documentation"))
print ("Modification of the code provided with exciting DFT code. The alternative definition of")
print ("deformation for cubic system has been implemented. The code works with exciting <input.xml>")
print ("but it can be converted to POSCAR format by <atomsk> code.")
print (CRED + " --> List of deformation codes for strains in Voigt notation, ONLY for Cubic. " + CEND )
print (CRED + " --> Values taken from the IRelast package paper and        "  + CEND )  
print (" --> Feng, W., Effects of short-range order on the magnetic and mechanical properties ")
print ("     of FeCoNi(AlSi)x high entropy alloys. Metals, 7(11), 482 (2017). ")
print (" --> https://doi.org/10.1016/j.jallcom.2017.10.139")
print ('{:_^80s}'.format("END of Documentation"))
print ("                        |η[0]    η[5]/2  η[4]/2| ")
print ("                    η = |η[5]/2  η[1]    η[3]/2| ")
print ("                        |η[4]/2  η[3]/2  η[2]  | ")  
print ("                               D' = I + η")       
print ("------------------------------------------------------------------------")		
print (CYEL + " 0 => (η, η,           η,0,0,0) | volumetric strain [3C11+6C12=9B0]" + CEND)
print (CYEL + " 1 => (η,-η,1/1-η**2 - 1,0,0,0) | volume-conserving orthorhombic distortion 2[C11-C12]" + CEND)
print (CYEL + " 2 => (0, 0,1/1-η**2 - 1,0,0,2η)| volume-conserving monoclinic distortion [2C44]" + CEND)
print ("------------------------------------------------------------------------")           
print ( "="*80 )

maximum_strain = float( input("Enter maximum Lagrangian strain [1-10%] >>>> ") )
if (1 < maximum_strain or maximum_strain < 0): 
    sys.exit("ERROR: Maximum Lagrangian strain is out of range [0-1]!\n")
strain_points = int( input("Enter # of strain values (odd preferably) >>>> ") )
tmp = int ( floor(strain_points/2) )
print("The deformation range is [-{},{}]".format(tmp, tmp) )
if (3 > strain_points or strain_points > 99): 
    sys.exit("ERROR: Number of strain values is out of range [3-99]!\n")
		
deformation_code = int(input("Enter deformation code >>>> "))
print ("------------------------------------------------------------------------")
if (0 > deformation_code or deformation_code > 3): 
    sys.exit("ERROR: Deformation code is out of range [0-3]!\n")

if (deformation_code == 0 ): dc='EEE000'
if (deformation_code == 1 ): dc='EEE000'
if (deformation_code == 2 ): dc='00E00E'

#-------------------------------------------------------------------------------

input_obj = open("input.xml","r")
input_doc = etree.parse(input_obj)
input_rut = input_doc.getroot()
 
xml_scale = input_doc.xpath('/input/structure/crystal/@scale')
if (xml_scale == []): 
    ref_scale = 1.0
else: 
    ref_scale = float(xml_scale[0])

str_stretch = input_doc.xpath('/input/structure/crystal/@stretch')
if (str_stretch ==[]): 
    xml_stretch = [1.,1.,1.]
else: 
	xml_stretch=np.array( float(str_stretch[i].split()) for i in range(str_stretch) )

lst_basevect = input_doc.xpath('//basevect/text()')
xml_basevect = []
for ind_basevect in lst_basevect:
	l = float(ind_basevect.split()[0]), float(ind_basevect.split()[1]), float(ind_basevect.split()[2])
	#l = np.array(l); print (l)
	xml_basevect.append( l )

axis_matrix = np.array(xml_basevect) 
determinant = np.linalg.det(axis_matrix)
volume = abs(determinant * ref_scale**3 * xml_stretch[0]*xml_stretch[1]*xml_stretch[2])

work_directory = 'workdir'
if (len(sys.argv) > 1): work_directory = sys.argv[1]
if (os.path.exists(work_directory)): shutil.rmtree(work_directory)
os.mkdir(work_directory)
os.chdir(work_directory)

output_info =open('INFO-elastic-constants',"w")
output_info.write("\n")
output_info.write("Maximum Lagrangian strain       = {}\n".format(maximum_strain))
output_info.write("Number of strain values         = {}\n".format(strain_points ))
output_info.write("Volume of equilibrium unit cell = {} [a.u]^3\n".format(volume))
output_info.write("Deformation code                = {}\n".format(deformation_code))
output_info.write("Deformation label               = {}\n".format(dc) )
output_info.close()

#-------------------------------------------------------------------------------

delta=strain_points-1
convert=1
if (strain_points <= 1):
    strain_points=1
    convert=-1
    delta=1
eta_step=2*maximum_strain/delta

#-------------------------------------------------------------------------------
t=1; tmp=-tmp; Vo=1
print ("{:12s} {:12.8s} {:14.8s} {:14.8s}".format("", "Vol_cell", "Vol_D'", "(V-Vo)/V" ))

for i in range(0,strain_points):
	eta=i*eta_step-maximum_strain*convert
	if (i+1 < 10): strainfile = 'strain-0'+str(i+1)
	else: strainfile = 'strain-'+str(i+1)
	
	output_str = open(strainfile,"w")
	output_str.write( "{:11.8f}\n".format(eta) )
	output_str.close()
	
	if (abs(eta) < 0): eta=0
	ep=eta
	if (eta < 0): em=abs(eta)
	else: em=-eta
#----------------------------------------------------------------------------
	e=[]
	for j in range(6):
		ev=0
		if  (dc[j:j+1] == 'E' ): ev=ep
		elif(dc[j:j+1] == 'e' ): ev=em
		elif(dc[j:j+1] == '0' ): ev=0
		else: print ("==> ", dc); sys.exit("ERROR: deformation code not allowed!") 
		e.append(ev) 
		#print (e)
#------------------------------------ DEFORMATION STRAIN MATRIX ---------------
	#     |e[0]    e[5]/2  e[4]/2|
	# D = |e[5]/2  e[1]    e[3]/2|
	#     |e[4]/2  e[3]/2  e[2]  |
	
	# BULK 9B = 3C11 + 6C12
	if (deformation_code == 0): 
		eta_matrix=np.matrix( [[ e[0], e[5], e[4]], [ e[5], e[1], e[3]],[ e[4], e[3], e[2]]] )
		
	# 2*[C11 - C12]
	if (deformation_code == 1): # 
		eta_matrix=np.matrix( [[ e[0], e[5], e[4]],[ e[5],-e[1], e[3]], [ e[4], e[3], 1/(1-e[2]**2) -1 ]] )
	
	# Here off diagonal terms are e/2 that is why we have 2*C44 otherwise we would have 4*C44
	if (deformation_code == 2): 
		eta_matrix=np.matrix( [[ e[0], e[5], e[4]],[ e[5], e[1], e[3]], [ e[4], e[3], 1/(1-e[2]**2) -1 ]] )
		
	one_matrix=np.identity(3)

#----------------------------------------------------------------------------
	norma=1.0; inorma=0 ; eps_matrix = eta_matrix
	if (linalg.norm(eta_matrix) > 0.7):sys.exit("ERROR: too large deformation!")
	
	while ( norma > 1.e-10 ):
		x=eta_matrix-0.5*dot(eps_matrix,eps_matrix)
		norma=linalg.norm(x-eps_matrix) 
		eps_matrix=x
		inorma=inorma+1			
	
	# deformation matrix D = 1 + eps = eps + 0.5*eps**2
	def_matrix = one_matrix + eps_matrix ; #print (def_matrix)	
	
	# transforming to new coordinates R' = (1 + eps).R
	new_axis_matrix=transpose(dot(def_matrix,transpose(axis_matrix)))
	nam=new_axis_matrix
	
	#                for debugging ONLY checking volume conserving condition		
###############################################################################################
	V = np.linalg.det(new_axis_matrix)
	V_def = np.linalg.det(def_matrix)
	
	# Defining Direction cosines:: Elastic anisotropy of crystals is important 
	# since it correlates with the possibility to induce micro-cracks in materials
	#n1 = np.linalg.norm(new_axis_matrix[0])
	#n2 = np.linalg.norm(new_axis_matrix[1])
	#n3 = np.linalg.norm(new_axis_matrix[2])
	#l1 = math.degrees(math.acos(np.dot(new_axis_matrix[0],np.transpose(new_axis_matrix[1]))/(n1*n2)))
	#l2 = math.degrees(math.acos(np.dot(new_axis_matrix[1],np.transpose(new_axis_matrix[2]))/(n2*n3)))
	#l3 = math.degrees(math.acos(np.dot(new_axis_matrix[2],np.transpose(new_axis_matrix[0]))/(n3*n1)))
	#print(l1,l2,l3)
	if (tmp == 0): # reference volume Vo
		Vo = abs(V)
		
	print ("{:2d}({:2d}) => {:10.6f} {:10.6f} {:14.6f}".format(t, tmp, abs(V), abs(V_def), (V-Vo)/V ) ) 
	#aa = np.linalg.det( (np.dot(axis_matrix,def_matrix) ) )
	#print(np.matmul( def_matrix,axis_matrix ) )
	#print ("{:8.6f}".format(aa) )
###############################################################################################	
#----------------------------------------------------------------------------
	
	xbv = input_doc.xpath('//crystal/basevect')
	fmt = '%22.16f'
	for j in range(3):
		xbv[j].text = str(fmt%nam[j,0])+str(fmt%nam[j,1])+str(fmt%nam[j,2])+" "
	if (i+1 < 10): outputfile = 'input-0'+str(i+1)+'.xml'
	else: outputfile = 'input-'+str(i+1)+'.xml'
	output_obj = open(outputfile,"wb")
	output_obj.write(etree.tostring(input_rut, method='xml',
																						pretty_print=True,
																						xml_declaration=False,
																						encoding='UTF-8'))
	output_obj.close()	
	t+=1; tmp +=1
#-------------------------------------------------------------------------------

os.chdir('../')
print("")


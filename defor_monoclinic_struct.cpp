#include <iostream>
#include <fstream>
#include <sstream>
#include <stdlib.h>
#include <iomanip>
#include <stdio.h>
#include <string>
#include <cmath>
#include <vector>
#include <algorithm>
#include <iterator>

using namespace std;

#define _USE_MATH_DEFINES

/*
**********************************DOCUMENTATION******************************
* Author: ASIF IQBAL BHATTI
* DATE : 24/10/2015
* VERSION: 1.0
* USAGE: icpc/icc/g++ -O3 defor.cpp (Intel/gcc compiler) tested.
* PURPOSE: Generating the deformed structures along each lattice vector
* ONLY for Monoclinic structure
* Reads only a POSCAR (VASP) file.
*****************************END OF DOCUMENTATION****************************
*/

template <class T>
bool from_string(T& t, const std::string& s, std::ios_base& (*f)(std::ios_base&))
{
  std::istringstream iss(s);
  return !(iss >> f >> t).fail();
}

int main()
{
	int array_size = 100000;
	int arr = 3;
    float a, b, l, ch, j, alpha, beta, gamma, fact, scale;
    float h1[3]={}, v1[3]={}, r1[3]={};
    float x, y, z, vol, cons, atnum;
    float latvec[arr][arr], xx;
    char  vectors[4] = "xyz";
    char  select;
    int loop = 0, loop2;
	string array[array_size];

    string adjust_vec;
    FILE * outputfile;
    std::cout.precision(12);
    std::ifstream indata;
    string filename;
    std::string line, coord;

    cout << "\n" << endl;
 	cout  << std::string(30, '-') <<  "\033[1;31m PROGRAM USAGE \033[0m" << std::string(30, '-') << endl;

    cout << endl;
    cout << "  [*] Deform volume keeping a:b:c ratio constant."  << endl;
    cout << "  [*] Vary lattice vector along each direction." << endl;
    cout << "  [*] Vary lattice while volume fixed." << "\n"  << endl;

 	std::cout << std::string(75, '-') << endl;
 	cout <<  "\033[1;32m ENTERING THE PROGRAM... \033[0m" << endl;

	cout << "  |  Enter the POSCAR filename : ";
	getline(cin, filename);
	indata.open(filename.c_str(), ios::in);

	if (indata.is_open()){
		cout << "  |  File opened succesfully: \n";
		cout << "  |  Reading ... ";
    	while (! indata.eof() )
			{
			getline (indata,line);
//			cout << line <<  endl;
			array[loop]=line;
			loop++;
			}
		cout <<"# of lines detected:   " << loop << endl;
    	std::cout << std::string(60, '-') << std::endl;

		indata.close();
		}
    else if(!indata.is_open())
		{
		cout << "   Unable to open file" << endl;
		exit(0);
		}
	cout << "The title of the file is:  \t" << array[0] <<endl;
	cout << "The scaling factor is:  \t"  << array[1] <<endl;

	for (loop2=2; loop2<5;loop2++)	{
		int i = loop2 - 1, j = 1;
		stringstream iss(array[loop2]);
		string buf;
		vector<string> tokens;
		while (iss >> buf)
			{
			tokens.push_back(buf);
  			if(from_string<float>(xx, std::string(buf), std::dec))
  				{
				latvec[i][j] = xx;
//				printf("%16.12f\t", latvec[i][j]* ::atof(array[1].c_str()));
// 				cout << i << "\t" << j << endl;
				j++;
  				}
			}
			std::cout << endl;
	}

//*************************** Display output*******************************
	for ( int k=1; k < 4; k++ )	{
		for ( int p=1; p < 4; p++ ) {
			printf("%16.12f\t",  latvec[k][p]* ::atof(array[1].c_str()));
			}
			std::cout << endl;
    	}

    cout << "\n"; cout << "Magnitude of latice vectors:   " << "\n" << endl;

	for ( int k=1; k < 4; k++ ) {
		for ( int p=1; p < 4; p++ ) {
			h1[k-1] = h1[k-1] + latvec[k][p]*latvec[k][p];
			}
	cout << "  |" << vectors[k-1]<< "| = " << "\t"; std::cout << sqrt(h1[k-1]) << endl;
	}

	for (int k = 0; k < 3; k++) { h1[k]=sqrt(h1[k]);}

//*****************************  Generating angles

	alpha = acos ((latvec[2][1]*latvec[3][1]+latvec[2][2]*latvec[3][2]+latvec[2][3]*latvec[3][3])/(h1[1]*h1[2])) *180/M_PI;
	printf("%4s %12.4f\n","  alpha:", alpha);
	beta  = acos ((latvec[1][1]*latvec[3][1]+latvec[1][2]*latvec[3][2]+latvec[1][3]*latvec[3][3])/(h1[0]*h1[2])) *180/M_PI;
	printf("%4s %12.4f\n","  beta :", beta);
	gamma = acos ((latvec[1][1]*latvec[2][1]+latvec[1][2]*latvec[2][2]+latvec[1][3]*latvec[2][3])/(h1[0]*h1[1])) *180/M_PI;
	printf("%4s %12.4f\n","  gamma:", gamma);
	printf("%4s %12.4f\n","  volume:", h1[0]*h1[1]*h1[2]*sin (beta*3.14/180) );
//
	cout << "\nSECTION A: To calculate lattice vectors with a given volume." << endl;
	cout << "This section concern with changing lattice constant if we were given volume." << endl;
	cout << "If volume --> delta(a, b, c)." << "\n" << endl;
	cout << "Enter yes or no to set lattice vectors:   "; cin>>adjust_vec;
	if (adjust_vec == "yes" or adjust_vec == "y" )
		{
		cout << "  Enter volume of a crystal to set lattice vectors:   " << "\t" << endl;cin>>vol;
//  Apply over here the scaling law.
// (ratio of volume)^3  = (ratio of lengths)^3
		fact = pow(vol/(h1[0]*h1[1]*h1[2]*sin (beta*3.14/180)),1.0/3.0);

		for (int k=1; k < 4; k++) {
			r1[k-1] = h1[k-1] * fact;
			printf("%12.4f\n", r1[k-1]);
		}
		printf("%4s %12.4f\n","  volume of the final structure:", r1[0]*r1[1]*r1[2]*sin (beta*3.14/180) );

		}
	else if (adjust_vec == "no" or adjust_vec == "n" )
		{ printf("Exiting from this section.");}

	cout << "Select an option along which to deform the crystal ... " << "\n" << endl;
	cout << "Type:   " << "x" << "   y " << "  z " << "  or  volume (v)" << endl << flush; cin >> select;
	cout << "Entered option is             :   " << "[" <<select << "]" << endl;
	cout << "Number of strain values (odd) :   "; cin>>b;
	if ( b <= 0) std::cout << "Strain values must be a positive integer" << "\n";
	cout << "Deformation value             :   "; std::cin >> a; cout <<"\n";
	if (a < 0) a = abs(a);
//****************** Check to see if the value is real or not***********************************
    while(std::cin.fail()) {
        std::cout << "Not real number:   " << std::endl;
        std::cin.clear();
        std::cin.ignore(256,'\n');
        std::cin >> a;
    	}

	ch = (2*a)/(b-1);
	j = 0;

	outputfile = fopen("deformationstruct.dat", "w");

	switch(select)
	{
		case 'x':
		case 'X':

		fprintf(outputfile, "%10s %14s %14s\n","X", "Y", "Z");
		fprintf(outputfile, "  ");; for (int i=0; i < 52; i++) {fprintf(outputfile,"_");}; fprintf(outputfile, "\n");
		for( int i = 0; i < b; i++ )
		{
			j = i * ch - a;
//			printf("%4d %10.4f %16.10f\n", i, j, j*h1[0] + h1[0]);
			fprintf(outputfile,"%16.10f %16.10f %16.10f\n", j*latvec[1][1] + latvec[1][1], j*latvec[1][2] + latvec[1][2], j*latvec[1][3] + latvec[1][3] );
			fprintf(outputfile,"%16.10f %16.10f %16.10f\n", latvec[2][1], latvec[2][2], latvec[2][3] );
			fprintf(outputfile,"%16.10f %16.10f %16.10f\n", latvec[3][1], latvec[3][2], latvec[3][3] );
			if (int (b/2.0) == i) cout << "    Original lattice vector from a file\n" << endl;
			fprintf(outputfile, "\n");
			printf("%16.10f %16.10f %16.10f\n", j*latvec[1][1] + latvec[1][1], j*latvec[1][2] + latvec[1][2], j*latvec[1][3] + latvec[1][3] );
			printf("%16.10f %16.10f %16.10f\n", latvec[2][1], latvec[2][2], latvec[2][3] );
			printf("%16.10f %16.10f %16.10f\n", latvec[3][1], latvec[3][2], latvec[3][3] );

			v1[0] = pow(j*latvec[1][1]  + latvec[1][1],2) + pow(j*latvec[1][2] + latvec[1][2],2) + pow(j*latvec[1][3] + latvec[1][3],2) ;
			v1[1] = pow(latvec[2][1],2) + pow(latvec[2][2],2) + pow(latvec[2][3],2);
			v1[2] = pow(latvec[3][1],2) + pow(latvec[3][2],2) + pow(latvec[3][3],2);
	        for (int k = 0; k < 3; k++) { v1[k]=sqrt(v1[k]);}
			printf("%6s %12.4f\n","  volume:", v1[0]*v1[1]*v1[2]*sin (beta*3.14/180) );
			printf("\n");
		}
		break;

		case 'y':
		case 'Y':

		fprintf(outputfile, "%10s %14s %14s\n","X", "Y", "Z");
		fprintf(outputfile, "  ");; for (int i=0; i < 52; i++) {fprintf(outputfile,"_");}; fprintf(outputfile, "\n");
		for( int i = 0; i < b; i++ )
		{
			j = i * ch - a;
//			printf("%4d %10.4f %16.10f\n", i, j, j*h1[1] + h1[1]);
			fprintf(outputfile,"%16.10f %16.10f %16.10f\n", latvec[1][1], latvec[1][2], latvec[1][3] );
			fprintf(outputfile,"%16.10f %16.10f %16.10f\n", j*latvec[2][1] + latvec[2][1], j*latvec[2][2] + latvec[2][2], j*latvec[2][3] + latvec[2][3] );
			fprintf(outputfile,"%16.10f %16.10f %16.10f\n", latvec[3][1], latvec[3][2], latvec[3][3] );
			if (int (b/2.0) == i) cout << "    Original lattice vector from a file\n" << endl;
			fprintf(outputfile, "\n");
			printf("%16.10f %16.10f %16.10f\n", latvec[1][1], latvec[1][2], latvec[1][3] );
			printf("%16.10f %16.10f %16.10f\n", j*latvec[2][1] + latvec[2][1], j*latvec[2][2] + latvec[2][2], j*latvec[2][3] + latvec[2][3] );
			printf("%16.10f %16.10f %16.10f\n", latvec[3][1], latvec[3][2], latvec[3][3] );

			v1[0] = pow(latvec[1][1],2) + pow(latvec[1][2],2) + pow(latvec[1][3],2);
			v1[1] = pow(j*latvec[2][1]  + latvec[2][1],2)     + pow(j*latvec[2][2] + latvec[2][2],2) + pow(j*latvec[2][3] + latvec[2][3],2);
			v1[2] = pow(latvec[3][1],2) + pow(latvec[3][2],2) + pow(latvec[3][3],2);
	        for (int k = 0; k < 3; k++) { v1[k]=sqrt(v1[k]);}
			printf("%6s %12.4f\n","  volume:", v1[0]*v1[1]*v1[2]*sin (beta*3.14/180) );
			printf("\n");
		}
		break;

		case 'z':
		case 'Z':

		fprintf(outputfile, "%10s %14s %14s\n","X", "Y", "Z");
		fprintf(outputfile, "  ");; for (int i=0; i < 52; i++) {fprintf(outputfile,"_");}; fprintf(outputfile, "\n");
		for( int i = 0; i < b; i++ )
		{
			j = i * ch - a;
//			printf("%4d %10.4f %16.10f\n", i, j, j*h1[2] + h1[2]);
			fprintf(outputfile,"%16.10f %16.10f %16.10f\n", latvec[1][1], latvec[1][2], latvec[1][3] );
			fprintf(outputfile,"%16.10f %16.10f %16.10f\n", latvec[2][1], latvec[2][2], latvec[2][3] );
			fprintf(outputfile,"%16.10f %16.10f %16.10f\n", j*latvec[3][1] + latvec[3][1], j*latvec[3][2] + latvec[3][2], j*latvec[3][3] + latvec[3][3] );
			if (int (b/2.0) == i) cout << "    Original lattice vector from a file\n" << endl;
			fprintf(outputfile, "\n");
			printf("%16.10f %16.10f %16.10f\n", latvec[1][1], latvec[1][2], latvec[1][3] );
			printf("%16.10f %16.10f %16.10f\n", latvec[2][1], latvec[2][2], latvec[2][3] );
			printf("%16.10f %16.10f %16.10f\n", j*latvec[3][1] + latvec[3][1], j*latvec[3][2] + latvec[3][2], j*latvec[3][3] + latvec[3][3] );

			v1[0] = pow(latvec[1][1],2) + pow(latvec[1][2],2) + pow(latvec[1][3],2);
			v1[1] = pow(latvec[2][1],2) + pow(latvec[2][2],2) + pow(latvec[2][3],2);
			v1[2] = pow(j*latvec[3][1] + latvec[3][1],2) + pow(j*latvec[3][2] + latvec[3][2],2) + pow(j*latvec[3][3] + latvec[3][3],2);
	        for (int k = 0; k < 3; k++) { v1[k]=sqrt(v1[k]);}
			printf("%6s %12.4f\n","  volume:", v1[0]*v1[1]*v1[2]*sin (beta*3.14/180) );
			printf("\n");
		}
		break;

		case 'v':
		case 'V':

		fprintf(outputfile, "%10s %14s %14s\n","X", "Y", "Z");
		fprintf(outputfile, "  ");; for (int i=0; i < 52; i++) {fprintf(outputfile,"_");}; fprintf(outputfile, "\n");
		for( int i = 0; i < b; i++ )
		{
			j = i * ch - a;
			fprintf(outputfile,"%16.10f %16.10f %16.10f\n", j*latvec[1][1] + latvec[1][1], j*latvec[1][2] + latvec[1][2], j*latvec[1][3] + latvec[1][3] );
			fprintf(outputfile,"%16.10f %16.10f %16.10f\n", j*latvec[2][1] + latvec[2][1], j*latvec[2][2] + latvec[2][2], j*latvec[2][3] + latvec[2][3] );
			fprintf(outputfile,"%16.10f %16.10f %16.10f\n", j*latvec[3][1] + latvec[3][1], j*latvec[3][2] + latvec[3][2], j*latvec[3][3] + latvec[3][3] );
			if (int (b/2.0) == i) cout << "    Original lattice vector from a file\n" << endl;
			fprintf(outputfile, "\n");
			printf("%16.10f %16.10f %16.10f\n",  j*latvec[1][1] + latvec[1][1], j*latvec[1][2] + latvec[1][2], j*latvec[1][3] + latvec[1][3]  );
			printf("%16.10f %16.10f %16.10f\n",  j*latvec[2][1] + latvec[2][1], j*latvec[2][2] + latvec[2][2], j*latvec[2][3] + latvec[2][3]  );
			printf("%16.10f %16.10f %16.10f\n",  j*latvec[3][1] + latvec[3][1], j*latvec[3][2] + latvec[3][2], j*latvec[3][3] + latvec[3][3]  );

			v1[0] = pow(j*latvec[1][1]  + latvec[1][1],2) + pow(j*latvec[1][2] + latvec[1][2],2) + pow(j*latvec[1][3] + latvec[1][3],2) ;
			v1[1] = pow(j*latvec[2][1]  + latvec[2][1],2) + pow(j*latvec[2][2] + latvec[2][2],2) + pow(j*latvec[2][3] + latvec[2][3],2);
			v1[2] = pow(j*latvec[3][1]  + latvec[3][1],2) + pow(j*latvec[3][2] + latvec[3][2],2) + pow(j*latvec[3][3] + latvec[3][3],2);
	        for (int k = 0; k < 3; k++) { v1[k]=sqrt(v1[k]);}
			printf("%8s %12.4f\n","  volume:", v1[0]*v1[1]*v1[2]*sin (beta*3.14/180) );
			printf("\n");
		}
		break;

		fclose(outputfile);

		default:
		cout << "Invalid option  " << endl;
	}
		return 0;


}

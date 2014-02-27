/*
	Fragbsse.cpp: Estimate BSSE (with estimate uncertainty) between two molecular fragments.
    Input is two xyzfiles, a QM method, and a confidence limit.
    Model parameters are found in bssedata.dat. Currently it contains parameters for:
        0: MP2/6-31G*
        1: MP2/aug-cc-pVDZ
    Confidence limit can be either 1 sigma (0) or 95% (1)

	e.g.  INPUT: ./fragbsse.exe 1.xyz 2.xyz 0 0 
	      OUTPUT: 2.0 0.5           
	Here we estimated the BSSE at MP2/6-31G* between 1.xyz and 2.xyz to be 2.0 kcal/mol with an error bar of 0.5 kcal/mol at 1 sigma confidence.

	Author:John Faver 2011
	See: Faver, J. C., Zheng, Z, Merz, K. M. Model for the Fast Estimation of Basis Set Superposition Error in Biomolecular Systems 
	     Journal of Chemical Physics, 2011, 135, 144110.10.1063/1.3641894

*/
#include "geom.cpp"
#include <iostream>
#include <string>
int main(int nargs, char* args[]){

	geom geom;
	int method;	
	int ci;

	if(nargs!=5){
		std::cout << "Usage: ./fragbsse.exe [*.xyz] [*.xyz] [method#] [confidence#]\n";
		std::cout << "METHODS\t0: MP2/6-31G*\tCONFIDENCE\t0: 68%\n\t1: MP2/aDZ\t\t\t1: 95% \n";
		exit(1);
	}

	std::string infilename=args[1];	
	molecule mol(infilename,"XYZ");
	std::string infilename2=args[2];
	molecule mol2(infilename2,"XYZ");
	method = atoi(args[3]);
	ci = atoi(args[4]);
	float* bsse = new float[2];	

	if (mol.charge < 0 or mol2.charge < 0)    			//anionic
		bsse = geom.getPscore(mol,mol2,method,0,ci);	
	else if(mol.charge > 0 or mol2.charge > 0)			//cationic
		bsse = geom.getPscore(mol,mol2,method,1,ci);	
	else if(geom.detectHbond(mol,mol2))				    //hbond
		bsse = geom.getPscore(mol,mol2,method,2,ci);	
	else	                							//vdw
		bsse = geom.getPscore(mol,mol2,method,3,ci);		

	std::cout <<bsse[0]<<" "<< bsse[1]<< "\n";
	
	return 0;
}

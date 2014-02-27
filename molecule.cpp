
#include <fstream>
#include <iostream>
#include "atom.cpp"
#include <string>
#include <stdlib.h>

class molecule{

	public:
		atom* atomList;
		std::string filename,filetype;	
		char molname[20];
		int natoms;
		int charge;
		molecule();
		molecule(std::string filename, std::string filetype);
		void printInfo();
		void readXYZfile(std::string filename);
};
molecule::molecule(){
	charge=0;	
	natoms=0;
	filename="DEFAULT";
	strcpy(molname,filename.c_str());
	filetype="DEFAULT";
}
molecule::molecule(std::string filename,std::string filetype){
	if(filetype.find("XYZ")!=std::string::npos or filetype.find("xyz")!=std::string::npos)		
		readXYZfile(filename);
	else{
		std::cout << "ERROR: Expected xyz files. \n";
		exit(1);	
	}
}
void molecule::readXYZfile(std::string filename){
	filetype="XYZ";
	std::string templine;	
	std::ifstream xyzfile(filename.c_str());
	getline(xyzfile,templine);
	natoms = atoi(templine.c_str());			
	getline(xyzfile,templine);
	if(sscanf(templine.c_str(),"%s %d",molname,&charge)<2){
		std::cout << "ERROR: xyz file needs name and charge on line 2.\n";
		exit(1);	
	}
	atomList = new atom[natoms];
	for(int i=0;i<natoms;i++){
		getline(xyzfile,templine);
		if(sscanf(templine.c_str(), "%c %e %e %e", &atomList[i].atsym, &atomList[i].x, &atomList[i].y, &atomList[i].z ) <4){
			std::cout << "ERROR: Trouble reading "<<filename<<" in line "<<i+3<<"\n";
			exit(1);
		}
	}			
	xyzfile.close();
}

void molecule::printInfo(){
	std::cout << natoms << " atoms. "<< charge << " charge. \n";
	for(int i=0;i<natoms;i++){
		atomList[i].printInfo();	
	}
}


/*
	geom.cpp: Does various molecular geometry related stuff.
*/

#include "molecule.cpp"
#include <cmath>
#include <cstring>
#include <fstream>

class geom{
	public:
		float PI;
		geom();
		float rad2deg(float r);
		float magnitude(float* v);
		float dotProduct(float* v1, float* v2);
		float getDistance(atom &a1, atom &a2);	
		float getDistancesq(atom &a1, atom &a2);
		float getAngle(atom &a1, atom &a2, atom &a3);
		float getPscore(molecule &m1, molecule &m2,float a,float b,float c);
		float* getPscore(molecule &m1, molecule &m2,float a,float b,float c,int n,float shalf,float pav,float den,float t);
		float* getPscore(molecule &m1,molecule &m2,int method,int inttype, int ci);
		int   findNearestAtom(molecule &m, int id);
		int   findNearestHeavyAtom(molecule &m, int id);
		bool  detectHbond(molecule &m1, molecule &m2);
};
geom::geom(){
	PI = 3.1415926535;
}
float geom::rad2deg(float r){
	return r*180.0/PI;
}
float geom::dotProduct(float* v1, float* v2){
	return v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2];
}
float geom::magnitude(float* v){
	return	sqrt(dotProduct(v,v));
}
float geom::getDistance(atom &a1,atom &a2){
	return sqrt(pow(a1.x-a2.x,2)+pow(a1.y-a2.y,2)+pow(a1.z-a2.z,2));
}
float geom::getDistancesq(atom &a1,atom &a2){
	return pow(a1.x-a2.x,2)+pow(a1.y-a2.y,2)+pow(a1.z-a2.z,2);
}
float geom::getAngle(atom &a1,atom &a2,atom &a3){
	float v21[3];  
	float v23[3];
	v21[0] = a1.x-a2.x;
	v21[1] = a1.y-a2.y;
	v21[2] = a1.z-a2.z; 	
	v23[0] = a3.x-a2.x;
	v23[1] = a3.y-a2.y;
	v23[2] = a3.z-a2.z;
	return rad2deg(acos(dotProduct(v21,v23)/magnitude(v21)/magnitude(v23)));
}
float geom::getPscore(molecule &m1,molecule &m2,float a,float b,float c){
	float pscore = 0.0;	
	float tdist = 0.0;
	for (int i=0;i<m1.natoms;i++){
		if (m1.atomList[i].atsym=='H'){
			continue;
		}
		for (int j=0;j<m2.natoms;j++){
			if(m2.atomList[j].atsym=='H'){
				continue;
			}
			tdist = getDistancesq(m1.atomList[i],m2.atomList[j]);
				pscore+=exp(-c*tdist);
		}
	}	
	pscore=a+b*pscore;
	return pscore;
}
float* geom::getPscore(molecule &m1,molecule &m2,float a,float b,float c,int n,float shalf,float pav,float den,float t){
	float* pscore = new float[2];	
	float tdist = 0.0;
	pscore[0]=0.0;
	pscore[1]=0.0;
	for (int i=0;i<m1.natoms;i++){
		if (m1.atomList[i].atsym=='H'){
			continue;
		}
		for (int j=0;j<m2.natoms;j++){
			if(m2.atomList[j].atsym=='H'){
				continue;
			}
			tdist = getDistancesq(m1.atomList[i],m2.atomList[j]);
				pscore[0]+=exp(-c*tdist);
		}
	}	
	pscore[1]=t*shalf*sqrt(1.0+1.0/n+pow((pscore[0]-pav),2)/den);
	pscore[0]=a+b*pscore[0];
	return pscore;
}
float* geom::getPscore(molecule &m1,molecule &m2,int method,int inttype,int ci){
	float* pscore = new float[2];	
	float tdist = 0.0;
	pscore[0]=0.0;
	pscore[1]=0.0;
	float a,b,c,t1s,t95,shalf,pav,den;	
	int n;
	int methodscan=-1;
	std::ifstream datafile("bssedata.dat");
	std::string templine;	
	bool methodfound=false;
	while(getline(datafile,templine)){
		sscanf(templine.c_str(),"%*c %d %*s",&methodscan);	
		if(methodscan==method){
			methodfound=true;
			break;
		}
	}
	if(!methodfound)
		std::cout << "ERROR: BSSE data not found, weird stuff is going to happen.\n";
	for (int i=0;i<=inttype;i++){
		getline(datafile,templine);
	}
	sscanf(templine.c_str(),"%e %e %e %d %e %e %e %e %e",&a,&b,&c,&n,&shalf,&pav,&den,&t1s,&t95);
	datafile.close();
	switch(ci){
		case 0:
			return getPscore(m1,m2,a,b,c,n,shalf,pav,den,t1s);
		case 1:
			return getPscore(m1,m2,a,b,c,n,shalf,pav,den,t95);
		default:
			std::cout << "ERROR: Incorrect option for t-value. \n";
	}
}
int geom::findNearestAtom(molecule &m, int id){
	float dmin = 1000.0;
	float td   = 0.0;
	unsigned int idmin  = id;
	for(int i=0;i<m.natoms;i++){
		if(i == id)
			continue;
		td = getDistancesq(m.atomList[i],m.atomList[id]);
		if(td < dmin){
			dmin = td;
			idmin = i;
		}
	}
	if(idmin != id)
		return idmin;
	else{
		std::cout << "Error finding nearest atom to ";
		m.atomList[id].printInfo();
		std::cout << " \n";
		exit(1);
	}
}
int geom::findNearestHeavyAtom(molecule &m, int id){
	float dmin = 1000.0;
	float td   = 0.0;
	unsigned int idmin  = id;
	for(int i=0;i<m.natoms;i++){
		if(i == id || m.atomList[i].atsym == 'H')
			continue;
		td = getDistancesq(m.atomList[i],m.atomList[id]);
		if(td < dmin){
			dmin = td;
			idmin = i;
		}
	}
	if(idmin != id)
		return idmin;
	else{
		std::cout << "Error finding nearest heavy atom to ";
		m.atomList[id].printInfo();
		std::cout << " \n";
		exit(1);
	}
}
bool geom::detectHbond(molecule &m1, molecule &m2){
//returns true or false for presence of hbond. Ignores C-H-X bonds.
// Define Hbond as D'-D-H-A-A'
// distance between h and acceptor
	float dmax = 3.0;
	float dmaxsq = dmax*dmax;
// phi angle = H-A-A' 
	float phimin = 65.0;
	float phimax = 180.0;
// theta angle = D-H-A 
	float thetamin = 90.0;
	float thetamax = 180.0;
// Reminder: These requirements are VERY lenient. 
	int donorId;
	int acceptorPrimeId;
	float theta;
	float phi;
// Check if m1 is the donor.
	for(int i=0;i<m1.natoms;i++){
		if (m1.atomList[i].atsym != 'H')
			continue;
		for(int j=0;j<m2.natoms;j++){
			if(m2.atomList[j].atsym == 'H' || m2.atomList[j].atsym == 'C')
				continue;
			if(getDistancesq(m1.atomList[i],m2.atomList[j]) < dmaxsq){
				donorId = findNearestHeavyAtom(m1,i);
				if(m1.atomList[donorId].atsym == 'C')
					continue;
				acceptorPrimeId = findNearestHeavyAtom(m2,j);
				theta = getAngle(m1.atomList[donorId],m1.atomList[i],m2.atomList[j]);
				phi = getAngle(m1.atomList[i],m2.atomList[j],m2.atomList[acceptorPrimeId]);
				if(theta > thetamin && theta < thetamax && phi > phimin && phi < phimax)
					return true;	
			}
			
		}
	}	
// Check if m2 is the donor.
	for(int i=0;i<m2.natoms;i++){
		if (m2.atomList[i].atsym != 'H')
			continue;
		for(int j=0;j<m1.natoms;j++){
			if(m1.atomList[j].atsym == 'H' || m1.atomList[j].atsym == 'C')
				continue;
			if(getDistancesq(m2.atomList[i],m1.atomList[j]) < dmaxsq){
				donorId = findNearestHeavyAtom(m2,i);
				if(m2.atomList[donorId].atsym == 'C')
					continue;
				acceptorPrimeId = findNearestHeavyAtom(m1,j);
				theta = getAngle(m2.atomList[donorId],m2.atomList[i],m1.atomList[j]);
				phi = getAngle(m2.atomList[i],m1.atomList[j],m1.atomList[acceptorPrimeId]);
				if(theta > thetamin && theta < thetamax && phi > phimin && phi < phimax)
					return true;	
			}
			
		}
	}	
	return false;
}

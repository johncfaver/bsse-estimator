
#include <iostream>
#include <string>
#include <string.h>
class atom{
	public:
		float x,y,z;
		char atsym;		
		atom();
		atom(char tatsym,float tx,float ty,float tz);
		void printInfo();
};
atom::atom(){
	x=0.;
	y=0.;	
	z=0.;
	atsym='X';
}
atom::atom(char tatsym,float tx,float ty,float tz){
	x=tx;
	y=ty;
	z=tz;
	atsym=tatsym;
}
void atom::printInfo(){
	std::cout << atsym << " " << x << " "<< y << " " << z;
	std::cout << "\n";
}


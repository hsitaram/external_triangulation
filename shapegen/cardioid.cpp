#include<iostream>
#include<fstream>
#include<cmath>
#define N 101
#define PI 3.14159265
int main()
{
	int np;
	double a;
	double b;
	double dth;
	double x,y;
	std::ofstream outfile("polypoints.inp");
	
	a=0.2; b=0.1; 
	np=N;

	dth=2*PI/double(np-1);

	outfile<<np-1<<"\n";

	for(int i=0;i<np-1;i++)
	{
		x=a*(2*cos(i*dth)-cos(2.2*i*dth))-0.5;
		y=a*(2*sin(i*dth)-sin(2.2*i*dth))-1.0;
		outfile<<x<<"\t"<<y<<"\n";	
	}


	outfile.close();

	return(0);
}

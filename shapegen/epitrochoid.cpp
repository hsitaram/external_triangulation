#include<iostream>
#include<fstream>
#include<cmath>
#define N 41
#define PI 3.14159265
int main()
{
	int np;
	double R;
	double r,d;
	double dth;
	double x,y;
	std::ofstream outfile("polypoints.inp");
	
	R=0.2; r=0.5*R; d=0.5*r; 
	np=N;

	dth=2*PI/double(np-1);

	outfile<<np-1<<"\n";

	for(int i=0;i<np-1;i++)
	{
		x=(R+r)*cos(i*dth)-d*cos((R+r)/r*i*dth)-0.5;
		y=(R+r)*sin(i*dth)-d*sin((R+r)/r*i*dth)-1.0;
		outfile<<x<<"\t"<<y<<"\n";	
	}


	outfile.close();

	return(0);
}

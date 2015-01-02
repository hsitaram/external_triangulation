#include<iostream>
#include<fstream>
#define N 64
int main()
{
	int np;
	double xmin,xmax;
	double ymin,ymax;
	double dx;
	std::ofstream outfile("polypoints.inp");

	xmin=0.0; xmax=1.0; ymin=0.0; ymax=1.0;

	dx=(xmax-xmin)/double(N-1);

	np=4*(N-1);

	outfile<<np<<"\n";

	for(int i=1;i<N;i++)
	{
		outfile<<xmin+dx*i<<"\t"<<ymin<<"\n";	
	}

	for(int i=1;i<N;i++)
	{
		outfile<<xmax<<"\t"<<ymin+dx*i<<"\n";
	}

	for(int i=1;i<N;i++)
	{
		outfile<<xmax-dx*i<<"\t"<<ymax<<"\n";
	}

	for(int i=1;i<N;i++)
	{
		outfile<<xmin<<"\t"<<ymax-dx*i<<"\n";
	}

	outfile.close();

	return(0);
}

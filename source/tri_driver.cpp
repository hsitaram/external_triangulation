#include"polygon.h"
#include"triangulation.h"
#include<cstdlib>

int main(int argc,char *argv[])
{
	int n,npoly;
	double *p;
	polygon *pgons;
	triangulation triobj;
	std::string polyfname;
	char num[20];
	std::string strnum;
	double minlength;
	double globalavglength;

	std::ifstream infile;
	infile.open("polypoints.inp");

	infile>>npoly;

	pgons = new polygon[npoly];
	globalavglength = 0.0;

	for(int i=0;i<npoly;i++)
	{
		infile>>n;
		p = new double[DIM*n];
		for(int j=0;j<n;j++)
		{
			infile>>p[DIM*j]>>p[DIM*j+1];
		}

		pgons[i].assignpolypoints(p,n);
		sprintf(num,"%3.3d",i);
		strnum=num;
		polyfname="poly_"+strnum+".dat";
		pgons[i].printpoints(polyfname.c_str());

		globalavglength += pgons[i].avglength;

		delete(p);
	}

	infile.close();
	globalavglength = globalavglength/double(npoly);

	triobj.setpolydomain(pgons,npoly);
	std::cout<<"set poly domain\n";
	triobj.generate_external_triangulation();
	triobj.centroidinsert(globalavglength);
	//triobj.bwalgorithm(pgon.minlength);
	//triobj.uniformly_distribute_points(30,30);

	//triobj.printtridata();
	triobj.printtrianglesgnuplot();
	triobj.printtrianglesvtu();

	return(0);
}

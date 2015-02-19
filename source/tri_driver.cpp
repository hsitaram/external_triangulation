#include"polygon.h"
#include"triangulation.h"
#include"laplace_solver.h"
#include<cstdlib>

int main(int argc,char *argv[])
{
	int n,npoly;
	double *p;
	polygon *pgons;
	std::string polyfname;
	char num[20];
	std::string strnum;
	double minlength;
	double globalavglength;

	triangulation *triobj;
	laplace_solver *solverobj;


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

	triobj=new triangulation;
	triobj->setpolydomain(pgons,npoly);
	std::cout<<"set poly domain\n";
	triobj->generate_external_triangulation();
	triobj->centroidinsert(globalavglength*2);
	triobj->printtrianglesgnuplot();
	triobj->printtridata();

	solverobj=new laplace_solver;
	solverobj->assignelements(triobj,2.0,0.0);
	solverobj->solve();
	//triobj.printtrianglesvtu();

	return(0);
}

#include"edge.h"
#define INFTY 3e8
#ifndef TRIANGLE_H
#define TRIANGLE_H
class triangle
{
	double x1,y1;
	double x2,y2;
	double x3,y3;

	int pointids[3];

	void findcircumcenter(std::string str="nil");
	void findarea();
	void findavgsidelength();

	public:
	
	double circ_x,circ_y,circrad;
	double area,avglen;
	int edgeids[3];
	double centroid[2];
	int edgepos;
	
	void setvertices(double px1,double py1,double px2,double py2,
			double px3,double py3,int pid1,int pid2,int pid3,std::string str="nil")
	{
		x1=px1; x2=px2; x3=px3; 
		y1=py1; y2=py2; y3=py3;
		pointids[0]=pid1; 
		pointids[1]=pid2; 
		pointids[2]=pid3;
		
		edgepos=0;
		findcircumcenter(str);
		findarea();
		findavgsidelength();
		centroid[0]=0.33333*(x1+x2+x3);
		centroid[1]=0.33333*(y1+y2+y3);
	}

	void getvertices(double *points)
	{
		points[0]=x1; points[1]=y1;
		points[2]=x2; points[3]=y2;
		points[4]=x3; points[5]=y3;
	}

	void getotheredgeids(int me,int other[2]);

	int * getpointids()
	{ return(pointids); }

	void getpointids(int tp[3])
	{ tp[0]=pointids[0]; tp[1]=pointids[1]; tp[2]=pointids[2]; }

	void getedges(edge &e1,edge &e2,edge &e3);
	int getoppositepoint(edge e);
	void printedges(std::ofstream &outfile);
	bool ispointinsidecircumcircle(double px,double py);
	bool ispointinsidetriangle(double px,double py);
	bool getotherpoints(int idin,int &idout1,int &idout2);
};
#endif

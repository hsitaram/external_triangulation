#include<iostream>
#include<cmath>
#include<vector>
#include<fstream>
#define TOL 0.00000001
#define DIM 2

#ifndef EDGE_H
#define EDGE_H
class edge
{

	private:

		double start[DIM];
		double end[DIM];
		double length;


	public:
		int adjcellids[2];
		int nodeids[2];
		double A,B,C;

		void setedgeparams(double x1,double y1,double x2,double y2,int id1=-1,int id2=-1)
		{ 
			start[0] = x1; start[1] = y1;
			end[0]   = x2; end[1]  =  y2;
			length = sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1));
			nodeids[0]=id1;    nodeids[1]=id2;
			adjcellids[0]=-1;  adjcellids[1]=-1;
			A = (y2-y1); B=-(x2-x1); C=y1*(x2-x1)-x1*(y2-y1);
		}

		void getendpoints(double *first,double *last)
		{
			first[0] = start[0]; first[1] = start[1];
			last[0]  = end[0];    last[1] = end[1];
		}
		void getendpoints(double &x1,double &y1,double &x2,double &y2)
		{
			x1=start[0]; y1=start[1]; x2=end[0]; y2=end[1];
		}

		double getlength(){return length;}

		void substitutelrcell(int oldcellnum,int newcellnum);
		bool edgeintersect(edge e);
		bool edgeintersect(edge e,double &intx,double &inty);
		bool isitsame(edge e);

		void printendpoints()
		{
			std::cout<<"edge end points:"<<start[0]<<"\t"<<
				start[1]<<"\t"<<end[0]<<"\t"<<end[1]<<"\n";
		}

};
#endif

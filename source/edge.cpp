#include"edge.h"

//==========================================================
bool edge::edgeintersect(edge e)
{
	double intx,inty;
	bool intersectflag;

	intersectflag=edgeintersect(e,intx,inty);
	return(intersectflag);
}
//==========================================================
bool edge::edgeintersect(edge e,double &intx,double &inty)
{
	double x1,y1;
	double x2,y2;
	double x3,y3;
	double x4,y4;

	double px,py;
	double dtr,ntr1,ntr2;
	double dotp1,dotp2;
	bool intersectflag;

	x1=start[0]; y1=start[1];
	x2=end[0];   y2=end[1];

	e.getendpoints(x3,y3,x4,y4);

	dtr = (x1-x2)*(y3-y4)-(y1-y2)*(x3-x4);

	ntr1 = (x1*y2-y1*x2)*(x3-x4) - (x1-x2)*(x3*y4-y3*x4);
	ntr2 = (x1*y2-y1*x2)*(y3-y4) - (y1-y2)*(x3*y4-y3*x4);

	if(fabs(dtr) < TOL)
	{
		if(fabs(ntr1) < TOL && fabs(ntr2) < TOL)
		{
			intersectflag=true;
		}	
		else
		{
			intersectflag=false;
		}
	}
	else
	{
		px = ntr1/dtr; py=ntr2/dtr;
		intx = px; inty = py;

		dotp1=(px-x1)*(px-x2)+(py-y1)*(py-y2);
		dotp2=(px-x3)*(px-x4)+(py-y3)*(py-y4);

		if(dotp1<=TOL && dotp2<=TOL)
		{
			intersectflag=true;
		}
		else
		{
			intersectflag=false;
		}
	}

	return(intersectflag);
	
}
//==========================================================
bool edge::isitsame(edge e)
{
	bool same;
	double x1,y1,x2,y2;
	double x3,y3,x4,y4;

	double midx1,midy1,midx2,midy2;
	double dtr;
	
	x1=start[0]; y1=start[1];
	x2=end[0];   y2=end[1];
	e.getendpoints(x3,y3,x4,y4);

	dtr = (x1-x2)*(y3-y4)-(y1-y2)*(x3-x4);

	midx1=0.5*(x1+x2); midx2=0.5*(x3+x4);
	midy1=0.5*(y1+y2); midy2=0.5*(y3+y4);

	same=false;

	if(fabs(dtr) < TOL)
	{
		if(fabs(midx1-midx2) < TOL && fabs(midy1-midy2) < TOL)
		{
			same=true;
		}	
	}

	return(same);
}
//==========================================================
void edge::substitutelrcell(int oldcellnum,int newcellnum)
{
	if(adjcellids[0] == oldcellnum)
	{
		adjcellids[0]=newcellnum;
	}
	else
	{
		adjcellids[1]=newcellnum;
	}
}
//==========================================================

#include"edge.h"

//===============================================
double findanglebtwedges(edge e1,edge e2)
{
	double e1x1,e1y1,e1x2,e1y2;
	double e2x1,e2y1,e2x2,e2y2;

	double vec1[2],vec2[2],dotp,vec1mag,vec2mag,angle;

	e1.getendpoints(e1x1,e1y1,e1x2,e1y2);	
	e2.getendpoints(e2x1,e2y1,e2x2,e2y2);

	vec1[0]=(e1x2-e1x1);
	vec1[1]=(e1y2-e1y1);
	vec2[0]=(e2x2-e2x1);
	vec2[1]=(e2y2-e2y1);

	vec1mag=sqrt(vec1[0]*vec1[0]+vec1[1]*vec1[1]);
	vec2mag=sqrt(vec2[0]*vec2[0]+vec2[1]*vec2[1]);

	dotp=vec1[0]*vec2[0]+vec1[1]*vec1[1];

	angle=acos(dotp/vec1mag/vec2mag);
	
	return(angle);
}
//===============================================
double findareaoftriangle(double x1,double y1,double x2,double y2,double x3,double y3)
{
 return(0.5*(x1*(y2-y3)+x2*(y3-y1)+x3*(y1-y2)));
}
//===============================================

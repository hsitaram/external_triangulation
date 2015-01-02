#include"triangle.h"

//==============================================================================================
void triangle::printedges(std::ofstream &outfile)
{
	outfile<<x1<<"\t"<<y1<<"\n"<<x2<<"\t"<<y2<<"\n\n";
	outfile<<x2<<"\t"<<y2<<"\n"<<x3<<"\t"<<y3<<"\n\n";
	outfile<<x3<<"\t"<<y3<<"\n"<<x1<<"\t"<<y1<<"\n\n";
}
//==============================================================================================
void triangle::findarea()
{
	area=fabs(0.5*(x1*(y2-y3)+x2*(y3-y1)+x3*(y1-y2)));
}
//==============================================================================================
void triangle::findavgsidelength()
{

	avglen  = sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2));
	avglen += sqrt((x2-x3)*(x2-x3)+(y2-y3)*(y2-y3));
	avglen += sqrt((x3-x1)*(x3-x1)+(y3-y1)*(y3-y1));

	avglen = 0.333333*avglen;
}
//==============================================================================================
void triangle::findcircumcenter(std::string str)
{

	double Sx,Sy,a;

	Sx = 0.5 *( (x1*x1+y1*y1)*(y2-y3) + (x2*x2+y2*y2)*(y3-y1) + (x3*x3+y3*y3)*(y1-y2) );
	Sy =-0.5 *( (x1*x1+y1*y1)*(x2-x3) + (x2*x2+y2*y2)*(x3-x1) + (x3*x3+y3*y3)*(x1-x2) );

	a = x1*(y2-y3)+x2*(y3-y1)+x3*(y1-y2);

	if(fabs(a)<=TOL)
	{
		std::cout<<str<<"\n";
		std::cout<<"ERROR:0 area for triangle\n";
		std::cout<<x1<<" "<<y1<<" "<<x2<<" "<<y2<<" "<<x3<<" "<<y3<<"\n";
		std::cout<<"a="<<a<<"\n";
		std::cout<<"SX="<<Sx<<"\n";
		std::cout<<"SY="<<Sy<<"\n";
		circ_x = INFTY; circ_y=INFTY;
		circrad = INFTY;
	}
	else
	{
		circ_x = Sx/a; circ_y=Sy/a;
		circrad = sqrt((x1-circ_x)*(x1-circ_x)+(y1-circ_y)*(y1-circ_y));	
	}
}
//==============================================================================================
void triangle::getedges(edge &e1,edge &e2,edge &e3)
{
	e1.setedgeparams(x1,y1,x2,y2,pointids[0],pointids[1]);	
	e2.setedgeparams(x2,y2,x3,y3,pointids[1],pointids[2]);	
	e3.setedgeparams(x3,y3,x1,y1,pointids[2],pointids[0]);

}
//==============================================================================================
bool triangle::ispointinsidecircumcircle(double px,double py)
{
	bool inside;
	double dist2;

	if(circrad == INFTY)
	{
		inside=true;
	}
	else
	{
		dist2 = (px-circ_x)*(px-circ_x)+(py-circ_y)*(py-circ_y);
		inside=(dist2<circrad*circrad)?true:false;
	}
	return(inside);
}
//==============================================================================================
int triangle::getoppositepoint(edge e)
{
	int id1,id2,i;

	id1=e.nodeids[0]; id2=e.nodeids[1];

	for(i=0;i<3;i++)
	{
		if(pointids[i]!=id1 && pointids[i]!=id2)
		{
			break;
		}
	}	

	return(pointids[i]);
}
//==============================================================================================
void triangle::getotheredgeids(int me,int other[2])
{
	int pos=0;

	for(int i=0;i<3;i++)
	{
		if(edgeids[i] != me)
		{
			if(pos==2)
			{
				std::cout<<"error:no matching edge in getotheredgeids\n";
				std::cout<<"me:"<<me<<"\n";
			}
			else
			{
				other[pos++]=edgeids[i];
			}
		}
	}
}
//==============================================================================================
bool triangle::ispointinsidetriangle(double px,double py)
{
	bool pointoutside;

	edge p_to_C;
	edge e1,e2,e3;

	p_to_C.setedgeparams(px,py,centroid[0],centroid[1]);
	getedges(e1,e2,e3);

	pointoutside = e1.edgeintersect(p_to_C) | 
		e2.edgeintersect(p_to_C) | e3.edgeintersect(p_to_C);

	return(!pointoutside);

}
//==============================================================================================
bool triangle::getotherpoints(int idin,int &idout1,int &idout2)
{
	int pos;
	bool found;

	found  = false;
	idout1 = -1;
	idout2 = -1;
	
	for(int i=0;i<3;i++)
	{
		if(pointids[i] == idin)
		{
			found=true;
			pos=i;
			break;
		}	
	}

	if(found)
	{
		if(pos==0)
		{
			idout1=pointids[1]; idout2=pointids[2];
		}
		else if(pos==1)
		{
			idout1=pointids[0]; idout2=pointids[2];
		
		}
		else
		{
			idout1=pointids[0]; idout2=pointids[1];
		}
		return(found);
	}
	else
	{
		return(found);
	}
}
//==============================================================================================

#include"edge.h"

int main()
{
	edge e1,e2;

	e1.setedgeparams(0.0,0.0,0.0,1.0);
	e2.setedgeparams(1.0,1.0,0.0000001,0.0000001);


	if(e1.edgeintersect(e2))
	{
		std::cout<<"e1 intersects e2\n";
	}
	else
	{
		std::cout<<"e1 does not intersect e2\n";
	}

	return(0);
}

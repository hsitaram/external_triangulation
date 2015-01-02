#include"polygon.h"

//==============================================================================================
void polygon::assignpolypoints(double *p,int n)
{
	int nextp,myp;
	numpoints=n;
	allpoints.resize(DIM*numpoints);
	minlength =  INFTY;
	maxlength = -INFTY;
	avglength = 0.0;

	cent_x = 0.0;
	cent_y = 0.0;

	for(int i=0;i<n;i++)
	{
		allpoints[DIM*i]  = p[DIM*i];
		allpoints[DIM*i+1]= p[DIM*i+1];	
	}

	polyedges.resize(numpoints);

	for(int i=0;i<numpoints;i++)
	{
		myp = i;
		nextp=(i+1)%numpoints;

		polyedges[i].setedgeparams(allpoints[DIM*myp],allpoints[DIM*myp+1],
				allpoints[DIM*nextp],allpoints[DIM*nextp+1]);

		if(polyedges[i].getlength() < minlength)
		{
			minlength=polyedges[i].getlength();
		}
		if(polyedges[i].getlength() > maxlength)
		{
			maxlength=polyedges[i].getlength();
		}

		avglength += polyedges[i].getlength();

		cent_x += allpoints[DIM*myp];
		cent_y += allpoints[DIM*myp+1];
	}

	avglength = avglength/double(numpoints);

	cent_x = cent_x/double(numpoints);
	cent_y = cent_y/double(numpoints);

	earpoints.resize(0);
	reflexpoints.resize(0);
}
//==============================================================================================
void polygon::printpoints(std::string fname)
{
	std::ofstream outfile;

	outfile.open(fname.c_str());

	for(int i=0;i<numpoints;i++)
	{
		outfile<<allpoints[DIM*i]<<"\t"<<allpoints[DIM*i+1]<<"\n";
	}
	outfile<<allpoints[0]<<"\t"<<allpoints[1]<<"\n";

	outfile.close();
}
//==============================================================================================
void polygon::findearandreflexpoints()
{
	int ahead,behind,me;
	int npoints;
	double area;

	earpoints.resize(0);
	reflexpoints.resize(0);
	npoints=polyearcutpoints.size();

	for(int i=0;i<npoints;i++)
	{
		ahead=(i+1)%npoints;
		behind=(i+npoints-1)%npoints;

		ahead=polyearcutpoints[ahead];
		behind=polyearcutpoints[behind];
		me = polyearcutpoints[i];

		area = findareaoftriangle(allpoints[DIM*behind],allpoints[DIM*behind+1],
				allpoints[DIM*me],allpoints[DIM*me+1],
				allpoints[DIM*ahead],allpoints[DIM*ahead+1]);

		if(area > TOL)
		{
			earpoints.push_back(i);
		}
		else
		{
			reflexpoints.push_back(i);
		}
	}	 
}
//==============================================================================================
void polygon::cutear()
{
	int count,rfp,rfp1,rfp2,rfp3;
	triangle temp;
	edge tempedge;
	int p1,p2,p3,npoints,epindex;
	polytriangles.resize(0);
	polyearcutpoints.resize(numpoints);
	bool edgeintersectionflag;

	for(int i=0;i<numpoints;i++)
	{
		polyearcutpoints[i]=i;	
	}

	npoints=polyearcutpoints.size();

	while(npoints>4)
	{
		findearandreflexpoints();

		edgeintersectionflag=true;
		epindex=0;

		while(edgeintersectionflag)
		{
			edgeintersectionflag=false;
			//for triangle
			p1=(earpoints[epindex]+npoints-1)%npoints;
			p2=earpoints[epindex];
			p3=(earpoints[epindex]+1)%npoints;

			p1=polyearcutpoints[p1];
			p2=polyearcutpoints[p2];
			p3=polyearcutpoints[p3];

			tempedge.setedgeparams(allpoints[DIM*p1],allpoints[DIM*p1+1],
					allpoints[DIM*p3],allpoints[DIM*p3+1]);

			count=0;
			for(unsigned int i=0;i<polyedges.size();i++)
			{
				if(polyedges[i].edgeintersect(tempedge))
				{
					count=count+1;
				}
			}
			if(count > 4)
			{
				edgeintersectionflag=true;
			}
			epindex += 1;
		}
	

		temp.setvertices(allpoints[DIM*p1],allpoints[DIM*p1+1],
				allpoints[DIM*p2],allpoints[DIM*p2+1],
				allpoints[DIM*p3],allpoints[DIM*p3+1],p1,p2,p3);

		polytriangles.push_back(temp);

		polyearcutpoints.erase(polyearcutpoints.begin()+earpoints[epindex-1]);	
		npoints=polyearcutpoints.size();		
	}


	findearandreflexpoints();

	if(reflexpoints.size() == 0)
	{
		p1=polyearcutpoints[0];
		p2=polyearcutpoints[1];
		p3=polyearcutpoints[2];
	
		temp.setvertices(allpoints[2*p1],allpoints[2*p1+1],
			allpoints[2*p2],allpoints[2*p2+1],
			allpoints[2*p3],allpoints[2*p3+1],p1,p2,p3);

		polytriangles.push_back(temp);
		
		p1=polyearcutpoints[0];
		p2=polyearcutpoints[2];
		p3=polyearcutpoints[3];
		
		temp.setvertices(allpoints[2*p1],allpoints[2*p1+1],
			allpoints[2*p2],allpoints[2*p2+1],
			allpoints[2*p3],allpoints[2*p3+1],p1,p2,p3);

		polytriangles.push_back(temp);
	
	}
	else
	{
		rfp  = reflexpoints[0];
		rfp1 = (rfp+1)%4;
		rfp2 = (rfp+2)%4;
		rfp3 = (rfp+3)%4;

		p1=polyearcutpoints[rfp];
		p2=polyearcutpoints[rfp2];
		p3=polyearcutpoints[rfp3];
		
		temp.setvertices(allpoints[2*p1],allpoints[2*p1+1],
			allpoints[2*p2],allpoints[2*p2+1],
			allpoints[2*p3],allpoints[2*p3+1],p1,p2,p3);

		polytriangles.push_back(temp);
		
		p1=polyearcutpoints[rfp];
		p2=polyearcutpoints[rfp2];
		p3=polyearcutpoints[rfp1];
		
		temp.setvertices(allpoints[2*p1],allpoints[2*p1+1],
			allpoints[2*p2],allpoints[2*p2+1],
			allpoints[2*p3],allpoints[2*p3+1],p1,p2,p3);

		polytriangles.push_back(temp);
		
	}

	//std::cout<<"total number of triangles:"<<polytriangles.size()<<"\n";
}
//=================================================================================================
void polygon::printtriangles()
{
	std::ofstream outfile;
	outfile.open("triangles.dat");

	for(unsigned int i=0;i<polytriangles.size();i++)
	{
		polytriangles[i].printedges(outfile);
	}

	outfile.close();
}
//=================================================================================================
bool polygon::ispointinside(double px,double py)
{
	int numintersects;
	edge e;
	double maxx;
	double factor;
	bool inside;
	double intx,inty;
	bool intersects,pointexists;
	std::vector<double> intersectedpoints;
	int nintersectedpoints;
	double dist2;
	double outsidep_x,outsidep_y;
	double length,l_dash;

	
	maxx=-INFTY;
	factor=3;
	
	intersectedpoints.resize(0);


	length = sqrt((cent_x-allpoints[0])*(cent_x-allpoints[0])
		 + (cent_y-allpoints[1])*(cent_y-allpoints[1]));

	l_dash = factor*length;

	outsidep_x = allpoints[0] + l_dash*(allpoints[0]-cent_x)/length;
	outsidep_y = allpoints[1] + l_dash*(allpoints[1]-cent_y)/length;


	std::cout<<"seed point:"<<outsidep_x<<"\t"<<outsidep_y<<"\n";
	e.setedgeparams(px,py,outsidep_x,outsidep_y);
	

	numintersects=0;
	for(int i=0;i<numpoints;i++)
	{
		intersects=polyedges[i].edgeintersect(e,intx,inty);

		if(intersects)
		{
			nintersectedpoints=intersectedpoints.size()/2;
			pointexists=false;
		
			for(int k=0;k<nintersectedpoints;k++)
			{
				dist2 = (intx-intersectedpoints[2*k])*(intx-intersectedpoints[2*k]);
				dist2 += (inty-intersectedpoints[2*k+1])*(inty-intersectedpoints[2*k+1]);
		
				if(dist2 < TOL)
				{
					pointexists=true;
				}
			}

			if(!pointexists)
			{
				numintersects+=1;
				intersectedpoints.push_back(intx);
				intersectedpoints.push_back(inty);
			}
		}
	}

	//std::cout<<"no: of intersects:"<<numintersects<<"\n";

	if(numintersects%2 == 1)
	{
		inside=true;
	}
	else
	{
		inside=false;
	}

	return(inside);
}
//=================================================================================================

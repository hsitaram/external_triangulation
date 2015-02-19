#include"trielement.h"

//===========================================================
void trielement::setelementattribs(double coord[6],int gindices[3])
{
	globalnodeindices[0]=gindices[0];	
	globalnodeindices[1]=gindices[1];	
	globalnodeindices[2]=gindices[2];	

	x1 = coord[0]; y1=coord[1];
	x2 = coord[2]; y2=coord[3];
	x3 = coord[4]; y3=coord[5];

	dxdzeta = x2-x1;
	dxdeta  = x3-x1;
	dydzeta = y2-y1;
	dydeta  = y3-y1;

	jacdet = (x2-x1)*(y3-y1)-(x3-x1)*(y2-y1);

	dzetadx =  (y3-y1)/jacdet;
	dzetady = -(x3-x1)/jacdet;
	detadx  = -(y2-y1)/jacdet;
	detady  =  (x2-x1)/jacdet;
}
//===========================================================
void trielement::computelocalAX(double X[3],double AX[3])
{
	double dpsi_i_dx,dpsi_j_dx;
	double dpsi_i_dy,dpsi_j_dy;

	double zeta,eta;

	//random values, the derivatives
	//dont depend on it anyway
	zeta = 0.25;
	eta  = 0.25;

	//loop over test function
	for(int i=0;i<3;i++)
	{
		AX[i] = 0.0;

		dpsi_i_dx =  basisfuncder(i,0,zeta,eta)*dzetadx;
		dpsi_i_dx += basisfuncder(i,1,zeta,eta)*detadx;
		
		dpsi_i_dy  = basisfuncder(i,0,zeta,eta)*dzetady;
		dpsi_i_dy += basisfuncder(i,1,zeta,eta)*detady;
		
		//loop over trial function
		for(int j=0;j<3;j++)
		{
			dpsi_j_dx  =  basisfuncder(j,0,zeta,eta)*dzetadx;
			dpsi_j_dx +=  basisfuncder(j,1,zeta,eta)*detadx;
			
			dpsi_j_dy  =  basisfuncder(j,0,zeta,eta)*dzetady;
			dpsi_j_dy +=  basisfuncder(j,1,zeta,eta)*detady;

			AX[i] += -(dpsi_i_dx*dpsi_j_dx
					+dpsi_i_dy*dpsi_j_dy)*
					fabs(jacdet)*0.5*X[j];
		}

	}	
}
//===========================================================
void trielement::computelocaldiag(double diagterm[3])
{
	double dpsi_i_dx;
	double dpsi_i_dy;
	double zeta,eta;

	//random values, the derivatives
	//dont depend on it anyway
	zeta = 0.25;
	eta  = 0.25;

	//loop over test function
	for(int i=0;i<3;i++)
	{
		dpsi_i_dx = basisfuncder(i,0,zeta,eta)*dzetadx;
		dpsi_i_dx += basisfuncder(i,1,zeta,eta)*detadx;
		
		dpsi_i_dy  = basisfuncder(i,0,zeta,eta)*dzetady;
		dpsi_i_dy += basisfuncder(i,1,zeta,eta)*detady;

		//the 0.5 is from area of right angled triangle in
		//in master element.
		diagterm[i] = -(dpsi_i_dx*dpsi_i_dx + dpsi_i_dy*dpsi_i_dy)
				*fabs(jacdet)*0.5;
	}	

}
//===========================================================
double trielement::basisfuncval(int i,double zeta,double eta)
{
	double value=0.0;

	if(i == 0)
	{
		value = 1.0-zeta-eta;
	}
	else
		if(i == 1)
		{
			value = zeta;
		}
		else
			if(i == 2)
			{
				value = eta;
			}
			else
			{
				std::cerr<<"Only three basis functions "
					<<"for first order element\n";
			}

	return(value);
}
//===========================================================
double trielement::basisfuncder(int i,int dernum,double zeta,double eta)
{
	double derval=0.0;

	if(i == 0)
	{
		derval=(dernum==0)?-1.0:-1.0;
	}
	else
		if(i == 1)
		{
			derval=(dernum==0)?1.0:0.0;
		}
		else
			if(i == 2)
			{
				derval=(dernum==0)?0.0:1.0;
			}
			else
			{
				std::cerr<<"Only three basis functions "
					<<"for first order element\n";
			}

	return(derval);
}
//===========================================================

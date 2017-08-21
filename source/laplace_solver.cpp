#include"laplace_solver.h"

//===========================================================
void laplace_solver::assignelements(triangulation *tri_ptr,double velx,double vely)
{
	double coord[6];
	int pids[3];
	triptr=tri_ptr;
	int m,n;

	vx = velx;
	vy = vely;
	nelements=triptr->trianglelist.size();
	nnodes=triptr->numnodes;
	boundarynodeids=triptr->boundarynodeids;

	elementlist.resize(nelements);

	for(int i=0;i<nelements;i++)
	{
		triptr->trianglelist[i].getvertices(coord);
		triptr->trianglelist[i].getpointids(pids);
		elementlist[i].setelementattribs(coord,pids);
	}
	
	m_kspdim = 5;       //dimension of Krylov subspace
	m_numit	   = 400;            //number of restart iterations
	m_matsize  = nnodes;

	bvec = new double[nnodes];
	solnvec = new double[nnodes];

	n = nnodes;
	m = m_kspdim;

	m_kspvectors =  new double[n*(m+1)];
	m_Hessbergmat = new double[(m+1)*m];

	std::cout<<"finished assign elements\n";	
}
//===========================================================
void laplace_solver::findAX(double *AX,double *X,int n)
{
	int gindices[3];
	double localX[3];
	double localAX[3];
	int id;

	for(int i=0;i<n;i++)
	{
		AX[i]=0.0;
	}

	for(int i=0;i<nelements;i++)
	{
		elementlist[i].getglobalindices(gindices);
		
		localX[0]=X[gindices[0]];
		localX[1]=X[gindices[1]];
		localX[2]=X[gindices[2]];

		elementlist[i].computelocalAX(localX,localAX);
		
		AX[gindices[0]] += localAX[0];
		AX[gindices[1]] += localAX[1];
		AX[gindices[2]] += localAX[2];
	}

	//boundries
	//-----------

	//rectangular boundaries
	for(unsigned int i=0;i<triptr->boundarynodeids.size();i++)
	{
		id = triptr->boundarynodeids[i];
		AX[id] = X[id];
	}

	//body boundaries
	for(unsigned int i=0;i<triptr->bodynodeids.size();i++)
	{
		id = triptr->bodynodeids[i];
		AX[id] = X[id];
	}

}
//===========================================================
void laplace_solver::solve()
{
	double *x0;
	double *ax;
	double *minvx;

	x0 = new double[nnodes];
	ax = new double[nnodes];
	minvx = new double[nnodes];

	for(int i=0;i<nnodes;i++)
	{
		x0[i]=0.0;
	}

	findbvec();
	m_performgmres(bvec,x0,solnvec);

	jacobiprecond(minvx,x0,nnodes);

	triptr->printtrianglesvtu(solnvec);
}
//===========================================================
void laplace_solver::findbvec()
{
	double x,y;
	int id;

	for(int i=0;i<nnodes;i++)
	{
		bvec[i]=0.0;
	}

	for(unsigned int i=0;i<triptr->boundarynodeids.size();i++)
	{
		id = triptr->boundarynodeids[i];

		x = triptr->nodes[DIM*id];
		y = triptr->nodes[DIM*id+1];

		bvec[id] = vx*(y-triptr->domain_mid_y) - vy*(x-triptr->domain_mid_x);
	}

	for(unsigned int i=0;i<triptr->bodynodeids.size();i++)
	{
		id=triptr->bodynodeids[i];
		bvec[id] = 0.0;
	}
}
//===========================================================
void laplace_solver::noprecond(double *MinvX,double *X,int n)
{
	for(int i=0;i<n;i++)
	{
		MinvX[i]=X[i];
	}
}
//===========================================================
void laplace_solver::jacobiprecond(double *MinvX,double *X,int n)
{
	int gindices[3];
	double localdiag[3];
	double *globaldiag;
	int id;

	//note: can optimize this allocation
	//by doing it only once and making globaldiag 
	//a private variable
	globaldiag=new double[n];

	for(int i=0;i<n;i++)
	{
		globaldiag[i]=0.0;
	}

	for(int i=0;i<nelements;i++)
	{
		elementlist[i].getglobalindices(gindices);

		elementlist[i].computelocaldiag(localdiag);
		
		globaldiag[gindices[0]] += localdiag[0];
		globaldiag[gindices[1]] += localdiag[1];
		globaldiag[gindices[2]] += localdiag[2];
	}

	//boundries
	//-----------

	//rectangular boundaries
	for(unsigned int i=0;i<triptr->boundarynodeids.size();i++)
	{
		id = triptr->boundarynodeids[i];
		globaldiag[id]=1.0;
	}

	//body boundaries
	for(unsigned int i=0;i<triptr->bodynodeids.size();i++)
	{
		id = triptr->bodynodeids[i];
		globaldiag[id]=1.0;
	}

	//std::cout<<"globaldiag:";
	for(int i=0;i<n;i++)
	{
		//std::cout<<globaldiag[i]<<"\t";
		MinvX[i]=X[i]/globaldiag[i];
	}
	//std::cout<<"\n\n";

}
//===========================================================
void laplace_solver::precond(double *MinvX,double *X,int n)
{
 	//noprecond(MinvX,X,n);
	jacobiprecond(MinvX,X,n);
}
//===========================================================
double laplace_solver::m_findnorm(double *v1,int n)
{
	double norm=0;

	for(int i=0;i<n;i++)
	{
		norm=norm+v1[i]*v1[i];
	}

	return(sqrt(norm));
}
//========================================================================
double laplace_solver::m_innerproduct(double *v1,double *v2,int n)
{
	double innerprod=0;

	for(int i=0;i<n;i++)
	{
		innerprod=innerprod+v1[i]*v2[i];
	}

	return(innerprod);
}
//========================================================================
void laplace_solver::m_getkspvector(double *v1,int vecnum)
{
	int n;
	n=m_matsize;

	for(int i=0;i<n;i++)
	{
		v1[i]=m_kspvectors[vecnum*n+i];
	}
}
//========================================================================
void laplace_solver::m_setkspvector(double *vec,int vecnum)
{
	int n;
	n=m_matsize;

	for(int i=0;i<n;i++)
	{
		m_kspvectors[vecnum*n+i]=vec[i];
	}
}
//========================================================================
void laplace_solver::m_addvectors(double *v1,double *v2,double *v12,
		int n,double a,double b)
{
	for(int i=0;i<n;i++)
	{
		v12[i]=a*v1[i]+b*v2[i];
	}
}
//========================================================================
void laplace_solver::m_copyvector(double *v1,double *v2,int n) //(dest,source,size)
{
	for(int i=0;i<n;i++)
	{
		v1[i]=v2[i];
	}
}
//========================================================================
bool laplace_solver::m_arnoldialgorithm(double *v1)
{
	int m,i,j,index,n;
	double *Avj,*vj,*vi;
	double *wj,*tempvec;
	double *MinvAvj;

	bool lucky; //when norm becomes 0, KSP 
	//is no longer linearly independent.
	//we would have got the best solution.

	m = m_kspdim;
        n = m_matsize;

	Avj     = new double[n]();
	MinvAvj = new double[n]();
	vj      = new double[n]();
	vi      = new double[n]();
	wj      = new double[n]();
	tempvec = new double[n]();

	m_copyvector(vj,v1,n); 
	m_setkspvector(vj,0);

	lucky=false;
	for(j=0;j<m;j++)
	{
		m_getkspvector(vj,j);
		findAX(Avj,vj,n);
		precond(MinvAvj,Avj,n);

		m_copyvector(Avj,MinvAvj,n);
		//Avj is now M^-1 A vj
		//remember we are solving M^-1 A X = M^-1 b
		
		for(i=0;i<=j;i++)
		{
			index=i*m+j;
		        m_getkspvector(vi,i);
			m_Hessbergmat[index]=m_innerproduct(Avj,vi,n);
		}

		m_copyvector(wj,Avj,n);

		for(i=0;i<=j;i++)
		{
			index=i*m+j;
			m_getkspvector(vi,i);
			m_addvectors(wj,vi,tempvec,n,1.0,-m_Hessbergmat[index]);
			m_copyvector(wj,tempvec,n);
		}
		

		m_Hessbergmat[(j+1)*m+j] = m_findnorm(wj,n);

		if(m_Hessbergmat[(j+1)*m+j] > 0.0)
		{
			for(int i=0;i<n;i++)
			{
				wj[i]=wj[i]/m_Hessbergmat[(j+1)*m+j];
			}
		}
		else
		{
			if(m_Hessbergmat[(j+1)*m+j] != m_Hessbergmat[(j+1)*m+j])
			{
				std::cout<<"Nan detected in arnoldi algorithm\n";
			}
			lucky=true;
			break;
		}		

		m_setkspvector(wj,j+1);

	}	

	return(lucky);
}
//========================================================================
void laplace_solver::m_leastsqminimize(double *y,double beta)
{
	int m;
	double *beta_e1;
	double c,s,h_up,h_down,dtr;
	double val1,val2;

	m=m_kspdim;
	beta_e1 = new double[m+1]();

	beta_e1[0] = beta;

	//convert H into QR
	for(int i=0;i<m;i++)
	{
		h_up   = m_Hessbergmat[i*m + i    ];
		h_down = m_Hessbergmat[(i+1)*m + i];

		dtr = sqrt(h_up*h_up + h_down*h_down);
		
		c=h_up/dtr; s=h_down/dtr;

		for(int j=0;j<m;j++)
		{
			h_up   = m_Hessbergmat[i*m+j];
			h_down = m_Hessbergmat[(i+1)*m+j];

			//perform rotations
			//ith row
			m_Hessbergmat[i*m + j    ] =  c*h_up+s*h_down;
			//(i+1)th row
			m_Hessbergmat[(i+1)*m + j] = -s*h_up+c*h_down;

		}

		val1 =  c*beta_e1[i] + s*beta_e1[i+1];
		val2 = -s*beta_e1[i] + c*beta_e1[i+1];

		beta_e1[i]=val1; beta_e1[i+1]=val2;
	
	}


	// ||Hm y - beta e1|| = || QR y - Q Q^T beta e1||
	// || Q ( Ry - Q^T beta e1) || = || Ry - Q^T beta e1||

	//solve least squares problem
	y[m-1] = beta_e1[m-1]/m_Hessbergmat[(m-1)*m+(m-1)];

	for(int i=m-2;i>=0;i--)
	{
		y[i]=beta_e1[i];

		for(int j=i+1;j<m;j++)
		{
			y[i]=y[i]-m_Hessbergmat[i*m+j]*y[j];
		}

		y[i]=y[i]/m_Hessbergmat[i*m+i];
	}

}
//========================================================================
void laplace_solver::m_performgmres(double *b,double *x0,double *x)
{
	int n,m,it;
	double beta;
	double *v1;
	double *tempvec;
	double *v,*y;
	bool arnoldistopped;
	double *r,*r0,*Ax0,*Ax;
	double *Minvr;
	double minresidual;
	
	n = m_matsize;
	m = m_kspdim;
	it=0;
	minresidual=1e-10;
	beta=INFTY;
	
	v1      = new double[n]();
	tempvec = new double[n]();
	v       = new double[n]();
	r	= new double[n]();
	r0	= new double[n]();
	Ax0	= new double[n]();
	Ax	= new double[n]();
	Minvr   = new double[n]();
	
	y       = new double[m]();

	//finding r0
	findAX(Ax0,x0,n);
	m_addvectors(b,Ax0,r0,n,1.0,-1.0);
	precond(Minvr,r0,n);
	m_copyvector(r0,Minvr,n);

	//initial residual is r0=M^-1(b-Ax0)
	//we are solving M^-1 A X = b
	m_copyvector(r,r0,n);
	m_copyvector(x,x0,n);

	while(it<m_numit && beta>minresidual)
	{
		it++;
		std::cout<<"restart iteration:"<<it<<"\t";
		beta = m_findnorm(r,n);

		for(int i=0;i<n;i++)
		{
			v1[i]=r[i]/beta;
		}

		arnoldistopped = m_arnoldialgorithm(v1);

		if(arnoldistopped)
		{
			std::cout<<"lucky condition\n";
			break;
		}

		m_leastsqminimize(y,beta);

		for(int i=0;i<m;i++)
		{
			m_getkspvector(v,i);
			m_addvectors(x,v,tempvec,n,1.0,y[i]);
			m_copyvector(x,tempvec,n);
		}	

		//finding new residual
		findAX(Ax,x,n);
		m_addvectors(b,Ax,r,n,1.0,-1.0);
		precond(Minvr,r,n);
		m_copyvector(r,Minvr,n);

		std::cout<<"norm of residual:"<<m_findnorm(r,n)<<"\n";
	}

}
//========================================================================

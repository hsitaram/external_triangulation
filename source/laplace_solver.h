#include"triangulation.h"
#include"trielement.h"

class laplace_solver
{
	private:
		double vx,vy;
		double *solnvec;	
		double *bvec;
		int nnodes;
		int nelements;
		triangulation *triptr;
		std::vector<int> boundarynodeids;
		std::vector<trielement>elementlist;

		//gmres variables and functions
		int m_kspdim;
		int m_numit;
		int m_matsize;
		double *m_kspvectors;
		double *m_Hessbergmat;
		
		double m_innerproduct(double *v1,double *v2,int n);
		double m_findnorm(double *v1,int n);
		
		void   m_getkspvector(double *v1,int vecnum);
		void   m_setkspvector(double *vec,int vecnum);
		void   m_addvectors(double *v1,double *v2,double *v12,
				int n,double a,double b);
		void   m_copyvector(double *v1,double *v2,int n); //(dest,source,size)
		bool   m_arnoldialgorithm(double *v1);
		void   m_leastsqminimize(double *y,double beta);
		void   m_performgmres(double *b,double *x0,double *x);
		void findAX(double *AX,double *X,int n);
		void noprecond(double *MinvX,double *X,int n);
		void jacobiprecond(double *MinvX,double *X,int n);
		void precond(double *MinvX,double *X,int n);
		void findbvec();

	public:
		void assignelements(triangulation *tri_ptr,double velx,double vely);
		void solve();
		
};

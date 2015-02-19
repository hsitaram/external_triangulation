#include"triangle.h"

class trielement
{
	private:
		double x1,y1;
		double x2,y2;
		double x3,y3;
		double jacdet;

		double dxdzeta,dxdeta;
		double dydzeta,dydeta;

		double dzetadx,dzetady;
		double detadx,detady;
		int globalnodeindices[3];

		double basisfuncval(int i,double zeta,double eta);
		double basisfuncder(int i,int dernum,double zeta,double eta);

	public:
		void setelementattribs(double coord[6],int gindices[3]);
		void computelocalAX(double X[3],double AX[3]);
		void computelocaldiag(double diagterm[3]);
		void getglobalindices(int gindices[3])
		{
			gindices[0]=globalnodeindices[0];
			gindices[1]=globalnodeindices[1];
			gindices[2]=globalnodeindices[2];
		}

};

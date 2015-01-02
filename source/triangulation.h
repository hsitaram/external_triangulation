#include"triangle.h"
#include"polygon.h"
#include<algorithm>

class triangulation
{
	int numnodes;
	int numpoly;
	double xmin_d,ymin_d,xmax_d,ymax_d;
	std::vector<double> nodes;
	std::vector<edge> edgelist;
	std::vector<triangle> trianglelist;
	polygon *polydomain;
	std::vector<int> boundarynodeids;
	
	bool performflips();
	bool isedgepresent(edge e,int &pos);
	void addedgetolist(edge e,int trinum);
	void reordertriedges();
	void addnewpoint(int tid,double x=0.0,double y=0.0,bool addatcentroid=true);
	void formedgelist();
	void deletetrifromlist(std::vector<int> deltrilist);

	void testformedges();	
	void assign_nodes(std::vector<double> points)
	{ 
		nodes=points; 
		numnodes=nodes.size()/2;
	}
	void assign_triangles(std::vector<triangle> alltriangles);
	void insert_points_to_triangulation(std::vector<double> points,int np);

	void delete_interior_triangles(polygon domain);

	polygon rectdomain;
	
	public:

        void setpolydomain(polygon *p,int npoly)
	{
		polydomain = new polygon[npoly];

		for(int i=0;i<npoly;i++)
		{
			polydomain[i]=p[i];
		}

		numpoly    = npoly;
		boundarynodeids.resize(0);
	}


	void centroidinsert(double minsidelen);
	void generate_external_triangulation();
	void laplace_smooth();
	void uniformly_distribute_points(int nx,int ny);
	
	void printtridata();
	void printtrianglesgnuplot();
	void printtrianglesvtu();

};
#include"triangle.h"
#include"common.h"

#ifndef POLYGON_H
#define POLYGON_H
class polygon
{
	private:
	
	int numpoints;
	std::vector<int> polyearcutpoints;

	std::vector<int> earpoints;
	std::vector<int> reflexpoints;

	std::vector<edge> polyedges;

	void findearandreflexpoints();

	double cent_x,cent_y;

	public:
	
	std::vector<double> allpoints;
	std::vector<triangle> polytriangles;
	double minlength;
	double maxlength;
	double avglength;

	void assignpolypoints(double *p,int n);
	void printpoints(std::string fname);
	void cutear();
	void printtriangles();
	bool ispointinside(double px,double py);

};
#endif

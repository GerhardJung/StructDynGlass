#include "struct_voronoi.h"
#include "defs.h"
#include "pbc.h"
#include "voro++_2d.hh"
using namespace voro;

void eval_struct_voronoi(){

    container_poly_2d con(-boxL/2.0,boxL/2.0,-boxL/2.0,boxL/2.0,10,10,true,true,16);

    // Add 1000 random points to the container
	for(int i=0;i<N;i++) {
		double sigma = 0.5;
		if (type_data[i] == 1) sigma=0.44;
		if (type_data[i] == 2) sigma=0.47;
		con.put(i,xyz_inherent_data[i][0],xyz_inherent_data[i][1],sigma);
	}
	
	// Sum the Voronoi cell areas and compare to the container area
	double carea=boxL*boxL,varea=con.sum_cell_areas();
	printf("Total container area    : %g\n"
	       "Total Voronoi cell area : %g\n"
	       "Difference              : %g\n",carea,varea,varea-carea);
    
    // Do a custom computation on the Voronoi cells, printing the IDs,
	// positions, and Voronoi cell areas to a file
	//con.print_custom("%i %x %y %a %E %c","particles_random.out");
	c_loop_all_2d vl(con);
	voronoicell_neighbor_2d c;
	int ij,q;double *pp;
	double cx,cy;
	vector<int> vi;
	if(vl.start()) do if(con.compute_cell(c,vl)) {
			ij=vl.ij;q=vl.q;pp=con.p[ij]+con.ps*q;
			c.centroid(cx,cy);
			c.neighbors(vi);
			cout << con.id[ij][q] << " "<< *pp << " " << pp[0] << " " << c.perimeter() << " " << c.area() << " " << sqrt(cx*cx+cy*cy) << " " << vi.size() <<  std::endl;
	} while(vl.inc());
}
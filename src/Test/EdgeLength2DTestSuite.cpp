/*----------------------------------------------------------------------------*/
/*
 * EdgeLength2DTestSuite.cpp
 *
 *  Created on: 19 janv. 2011
 *      Author: ledouxf
 */
/*----------------------------------------------------------------------------*/
#include "EdgeLength2DTestSuite.h"
#include "Triton2/Core/Delaunay2D.h"
//#include "Triton2/Core/GeomToolKit.h"
//#include "Triton2/Core/RobustGeomToolKit.h"
#include "GMDS/IG/IGMesh.h"
#include "GMDSCEA/LimaWriter.h"
#include "GMDS/Utils/CommonTypes.h"
//#include "GMDSMeshTools/LaplacianSmoothing.h"
//#include "GMDSMeshTools/BoundaryOperator.h"
/*----------------------------------------------------------------------------*/
#include<vector>
#include <math.h>
/*----------------------------------------------------------------------------*/
using namespace triton;
using namespace gmds;
/*----------------------------------------------------------------------------*/
void EdgeLength2DTestSuite::testDoubleCircles(){
	std::vector<std::pair<math::Point,
				math::Point > > segments;
	std::vector<math::Point > pntsIN;
	std::vector<math::Point > pntsOUT;


	//INITIAL BOUNDARY CREATION BEGIN
	double pi = M_PI;
	for(double dtheta=0; dtheta<2.0*pi-10e-11; dtheta=dtheta+pi/32){
		double dx = cos(dtheta);
		double dy = sin(dtheta);
		pntsIN.push_back(math::Point(dx, dy));
		pntsOUT.push_back(math::Point(2*dx, 2*dy));
	}
	for(unsigned int i=0;i<pntsIN.size()-1;i++){
		segments.push_back(
				std::pair<math::Point,
				math::Point >(pntsIN[i],pntsIN[i+1]) );
	}
	segments.push_back(
			std::pair<math::Point,
			math::Point >(pntsIN[pntsIN.size()-1],pntsIN[0]) );

	for(unsigned int i=0;i<pntsOUT.size()-1;i++){
		segments.push_back(
				std::pair<math::Point,
				math::Point >(pntsOUT[i],pntsOUT[i+1]) );
	}
	segments.push_back(
			std::pair<math::Point,
			math::Point >(pntsOUT[pntsOUT.size()-1],pntsOUT[0]) );

	//INITIAL BOUNDARY CREATION END

	Delaunay2D<gmds::TCoord,GeomToolKit > ::OrientedBoundary bndIN(pntsIN);
	Delaunay2D<gmds::TCoord,GeomToolKit > ::OrientedBoundary bndOUT(pntsOUT);
	Delaunay2D<gmds::TCoord,GeomToolKit >  algo(bndOUT,&bndIN,1);
	algo.mesh(Delaunay2D<gmds::TCoord,GeomToolKit >::EDGE_LENGTH_STRATEGY);
	std::cout<<"Meshing done\n";

	IGMesh& m = algo.getMesh();
	std::cout<<"Nb faces: "<<m.getNbFaces()<<"\n";


//	int boundary = m.getNewMark();
//	gmds::BoundaryOperator<DIM2|N|F|F2N|F2F> bound_op(m);
//	bound_op.markBoundaryNodes(boundary);
//
//	gmds::LaplacianSmoothing<DIM2|N|F|F2N|F2F> smooth_op(m);
//	smooth_op.perform(100,boundary);
//	std::cout<<"smoothing done"<<std::endl;
//	m.unmarkAll(boundary);
//	m.freeMark(boundary);
//	std::cout<<"Nb faces: "<<m.getNbFaces()<<"\n";

	gmds::LimaWriter<IGMesh> writer(m);
	writer.write("delaunay_double_circles_EL",N|F);

	std::cout<<"VTK Writing done\n";
	CPPUNIT_ASSERT(1==1);
}
/*----------------------------------------------------------------------------*/

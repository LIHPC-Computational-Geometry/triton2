/*----------------------------------------------------------------------------*/
/*
 * Delaunay2DTestSuite.cpp
 *
 *  Created on: 19 janv. 2011
 *      Author: ledouxf
 */
/*----------------------------------------------------------------------------*/
#include "Delaunay2DTestSuite.h"
#include "Triton2/Core/Delaunay2D.h"
//#include "Triton/Core/GeomToolKit.h"
//#include "Triton/Core/RobustGeomToolKit.h"
#include "GMDS/IG/IGMesh.h"

#include "GMDS/Math/Point.h"
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
void Delaunay2DTestSuite::testInitDelaunay(){
	std::vector<gmds::math::Point> pnts;
	pnts.push_back(math::Point(1,1));
	pnts.push_back(math::Point(1,2));
	pnts.push_back(math::Point(2,1));
	pnts.push_back(math::Point(2,2));
	pnts.push_back(math::Point(1.7,1.3));
	pnts.push_back(math::Point(1.2,1.2));
	pnts.push_back(math::Point(2,1.4));
	pnts.push_back(math::Point(2,1.8));
	pnts.push_back(math::Point(2,1.6));

	Delaunay2D<gmds::TCoord,GeomToolKit > algo(pnts);
	algo.triangulate();

	IGMesh& m = algo.getMesh();
	LimaWriter<IGMesh> writer(m);
	writer.write("delaunay",N|F);
	CPPUNIT_ASSERT(1==1);
}
/*----------------------------------------------------------------------------*/
void Delaunay2DTestSuite::testDelaunaySquare(){
	std::vector<gmds::math::Point> pnts;

	for(int dx=0; dx<10; dx++)
		for(int dy= 0;dy<10;dy++)
				pnts.push_back(math::Point(dx, dy));

	Delaunay2D<gmds::TCoord,GeomToolKit > algo(pnts);
	algo.triangulate();


	IGMesh& m = algo.getMesh();
	LimaWriter<IGMesh> writer(m);
	writer.write("delaunay_square",N|F);
	CPPUNIT_ASSERT(1==1);
}
/*----------------------------------------------------------------------------*/
void Delaunay2DTestSuite::testDelaunayUnitSquare(){
	std::vector<gmds::math::Point> pnts;

//	for(int dx=0; dx<4; dx++)
//		for(int dy= 0;dy<4;dy++)
//				pnts.push_back(math::Point(dx, dy));
	for(int dx=0; dx<10; dx++)
			for(int dy= 0;dy<10;dy++)
				if(dx==0 || dx==9 || dy==0 || dy==9)
					pnts.push_back(math::Point(dx, dy));
	Delaunay2D<gmds::TCoord,GeomToolKit > algo(pnts);
	algo.meshNew();


	IGMesh& m = algo.getMesh();
	LimaWriter<IGMesh> writer(m);
	writer.write("delaunay_unit_square",N|F);
	CPPUNIT_ASSERT(1==1);
}
/*----------------------------------------------------------------------------*/
void Delaunay2DTestSuite::testDelaunayEmptySquare(){
	std::vector<gmds::math::Point> pnts;

	std::cout<<"TEST"<<std::endl;
	for(int dx=0; dx<10; dx++)
		for(int dy= 0;dy<10;dy++)
			if(dx==0 || dx==9 || dy==0 || dy==9)
				pnts.push_back(math::Point(dx, dy));

	Delaunay2D<gmds::TCoord,GeomToolKit >  algo(pnts);
	algo.mesh();


	IGMesh& m = algo.getMesh();
	LimaWriter<IGMesh> writer(m);
	writer.write("delaunay_empty_square",N|F);
	CPPUNIT_ASSERT(1==1);
}
/*----------------------------------------------------------------------------*/
void Delaunay2DTestSuite::testDelaunayCircle(){
	std::vector<gmds::math::Point> pnts;

	double pi = M_PI;
	for(double dtheta=0; dtheta<2.0*pi-10e-11; dtheta=dtheta+pi/32){
		double dx = cos(dtheta);
		double dy = sin(dtheta);
		pnts.push_back(math::Point(dx, dy));
	}


	Delaunay2D<gmds::TCoord,GeomToolKit > algo(pnts);
	algo.triangulate();

	IGMesh& m = algo.getMesh();
	LimaWriter<IGMesh> writer(m);
	writer.write("delaunay_circle",N|F);
	CPPUNIT_ASSERT(1==1);
}
/*----------------------------------------------------------------------------*/
void Delaunay2DTestSuite::testNonConvex(){
	std::vector<gmds::math::Point> pnts;

	int dx, dy;
	dy=0;
	for(dx=0; dx<10; dx++)
		pnts.push_back(math::Point(dx, dy));

	dy = 4;
	for(dx=0; dx<10; dx++)
		pnts.push_back(math::Point(dx, dy));

	dx=0;
	for(dy=1; dy<4; dy++)
		pnts.push_back(math::Point(dx, dy));

	dy = 1;
	for(dx=2; dx<10; dx++)
		pnts.push_back(math::Point(dx, dy));

	dy = 3;
	for(dx=2; dx<10; dx++)
		pnts.push_back(math::Point(dx, dy));

	pnts.push_back(math::Point(2,2));


	Delaunay2D<gmds::TCoord,GeomToolKit >  algo(pnts);
	algo.triangulate();

	IGMesh& m = algo.getMesh();
	LimaWriter<IGMesh> writer(m);
	writer.write("delaunay_non_convex",N|F);
	CPPUNIT_ASSERT(1==1);
}
/*----------------------------------------------------------------------------*/
void Delaunay2DTestSuite::testNonConvex2(){
	std::vector<gmds::math::Point> pnts;

	pnts.push_back(math::Point(0,0));
	pnts.push_back(math::Point(9,0));
	pnts.push_back(math::Point(9,1));
	pnts.push_back(math::Point(9,3));
	pnts.push_back(math::Point(9,4));
	pnts.push_back(math::Point(8,4));
	pnts.push_back(math::Point(5,4));
	pnts.push_back(math::Point(5.1,4));
	pnts.push_back(math::Point(3,4));
	pnts.push_back(math::Point(7.5,4));

	pnts.push_back(math::Point(2,1));
	pnts.push_back(math::Point(2,3));

	pnts.push_back(math::Point(2.1,3));
	pnts.push_back(math::Point(2.2,3));
	pnts.push_back(math::Point(2.5,3));


	pnts.push_back(math::Point(0,4));



	Delaunay2D<gmds::TCoord,GeomToolKit >  algo(pnts);
	algo.triangulate();

	IGMesh& m = algo.getMesh();
	LimaWriter<IGMesh> writer(m);
	writer.write("delaunay_non_convex2",N|F);
	CPPUNIT_ASSERT(1==1);
}
/*----------------------------------------------------------------------------*/
void Delaunay2DTestSuite::testNonConvex3(){
	std::vector<std::pair<math::Point,
	gmds::math::Point> > segments;

	gmds::math::Point P0 (0,0);
	gmds::math::Point P1 (9,0);
	gmds::math::Point P2 (9,1);
	gmds::math::Point P8 (9,3);
	gmds::math::Point P9 (9,4);
	gmds::math::Point P10(8,4);
	gmds::math::Point P13(5,4);
	gmds::math::Point P12(5.1,4);
	gmds::math::Point P14(3,4);
	gmds::math::Point P11(7.5,4);

	gmds::math::Point P3(2,1);
	gmds::math::Point P4(2,3);

	gmds::math::Point P5(2.1,3);
	gmds::math::Point P6(2.2,3);
	gmds::math::Point P7(2.5,3);


	gmds::math::Point P15(0,4);

	std::vector<gmds::math::Point> points;
	points.push_back(P0);
	points.push_back(P1);
	points.push_back(P2);
	points.push_back(P3);
	points.push_back(P4);
	points.push_back(P5);
	points.push_back(P6);
	points.push_back(P7);
	points.push_back(P8);
	points.push_back(P9);
	points.push_back(P10);
	points.push_back(P11);
	points.push_back(P12);
	points.push_back(P13);
	points.push_back(P14);
	points.push_back(P15);

//	Delaunay2D<gmds::TCoord,GeomToolKit >::OrientedBoundary boundary(points);
//	Delaunay2D<gmds::TCoord,GeomToolKit >  algo(boundary);
//	algo.triangulateConstrained();
//
//	gmds::Mesh<DIM2|N|F|F2N|F2F>& m = algo.getMesh();
//	gmds::VTKWriter<DIM2|N|F|F2N|F2F> writer(m);
//	writer.write("delaunay_non_convex3",N|F);
//	CPPUNIT_ASSERT(1==1);
}
/*----------------------------------------------------------------------------*/
void Delaunay2DTestSuite::testNonConvex4(){
	std::vector<gmds::math::Point> pnt;

	pnt.push_back(math::Point(0,0));
	pnt.push_back(math::Point(1,0));
	pnt.push_back(math::Point(2,0));
	pnt.push_back(math::Point(3,0));
	pnt.push_back(math::Point(4,0));
	pnt.push_back(math::Point(5,0));
	pnt.push_back(math::Point(6,0));
	pnt.push_back(math::Point(7,0));
	pnt.push_back(math::Point(8,0));
	pnt.push_back(math::Point(8.5,0));
	pnt.push_back(math::Point(9,0));
	pnt.push_back(math::Point(9,0.25));
	pnt.push_back(math::Point(9,0.5));
	pnt.push_back(math::Point(9,0.75));
	pnt.push_back(math::Point(9,1));
	pnt.push_back(math::Point(8.5,1));
	pnt.push_back(math::Point(8,1));
	pnt.push_back(math::Point(7,1));
	pnt.push_back(math::Point(6,1));
	pnt.push_back(math::Point(5,1));
	pnt.push_back(math::Point(4,1));
	pnt.push_back(math::Point(3,1));
	pnt.push_back(math::Point(2,1));
	pnt.push_back(math::Point(2,2));
	pnt.push_back(math::Point(2,3));
	pnt.push_back(math::Point(3,3));
	pnt.push_back(math::Point(4,3));
	pnt.push_back(math::Point(5,3));
	pnt.push_back(math::Point(6,3));
	pnt.push_back(math::Point(7,3));
	pnt.push_back(math::Point(8,3));
	pnt.push_back(math::Point(9,3));
	pnt.push_back(math::Point(9,3.25));
	pnt.push_back(math::Point(9,3.5));
	pnt.push_back(math::Point(9,3.75));
	pnt.push_back(math::Point(9,4));
	pnt.push_back(math::Point(8,4));
	pnt.push_back(math::Point(7,4));
	pnt.push_back(math::Point(6,4));
	pnt.push_back(math::Point(5,4));
	pnt.push_back(math::Point(4,4));
	pnt.push_back(math::Point(3,4));
	pnt.push_back(math::Point(2,4));
	pnt.push_back(math::Point(1,4));
	pnt.push_back(math::Point(0,4));
	pnt.push_back(math::Point(0,3));
	pnt.push_back(math::Point(0,2));
	pnt.push_back(math::Point(0,1));

	Delaunay2D<gmds::TCoord,GeomToolKit > ::OrientedBoundary boundary(pnt);
	Delaunay2D<gmds::TCoord,GeomToolKit >  algo(boundary);
	algo.mesh();
	IGMesh& m = algo.getMesh();
	LimaWriter<IGMesh> writer(m);
	writer.write("delaunay_non_convex4",N|F);
	CPPUNIT_ASSERT(1==1);
}
/*----------------------------------------------------------------------------*/
void Delaunay2DTestSuite::testDoubleCircles(){
	std::vector<std::pair<math::Point,
	gmds::math::Point> > segments;
	std::vector<gmds::math::Point> pntsIN;
	std::vector<gmds::math::Point> pntsOUT;

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
				gmds::math::Point>(pntsIN[i],pntsIN[i+1]) );
	}
	segments.push_back(
			std::pair<math::Point,
			gmds::math::Point>(pntsIN[pntsIN.size()-1],pntsIN[0]) );

	for(unsigned int i=0;i<pntsOUT.size()-1;i++){
		segments.push_back(
				std::pair<math::Point,
				gmds::math::Point>(pntsOUT[i],pntsOUT[i+1]) );
	}
	segments.push_back(
			std::pair<math::Point,
			gmds::math::Point>(pntsOUT[pntsOUT.size()-1],pntsOUT[0]) );


	Delaunay2D<gmds::TCoord,GeomToolKit > ::OrientedBoundary bndIN(pntsIN);
	Delaunay2D<gmds::TCoord,GeomToolKit > ::OrientedBoundary bndOUT(pntsOUT);
	Delaunay2D<gmds::TCoord,GeomToolKit >  algo(bndOUT,&bndIN,1);
//	algo.mesh();
//
//	gmds::Mesh<DIM2|N|F|F2N|F2F>& m = algo.getMesh();
//
//	int boundary_mark = m.getNewMark();
//	gmds::BoundaryOperator<DIM2|N|F|F2N|F2F> boundary_op(m);
//	boundary_op.markBoundaryNodes(boundary_mark);
//	gmds::Timer t3;
//	std::cout<<"Boundary marking time: "<<t3-t2<<"\n";
//
//	gmds::LaplacianSmoothing<DIM2|N|F|F2N|F2F> smoother(m);
//	smoother.perform(100, boundary_mark);
//	gmds::Timer t4;
//	std::cout<<"Laplacian smoothing time: "<<t4-t3<<"\n";
//
//	gmds::VTKWriter<DIM2|N|F|F2N|F2F> writer(m);
//	writer.write("delaunay_double_circles",N|F);
//	std::cout<<"VTK Writing time: "<<t5-t4<<"\n";
	CPPUNIT_ASSERT(1==1);
}
/*----------------------------------------------------------------------------*/
void Delaunay2DTestSuite::testSwapConvex(){
	std::vector<gmds::math::Point> pnts;

	Delaunay2D<gmds::TCoord,GeomToolKit >  algo(pnts);
	IGMesh& m = algo.getMesh();

	Node n0 = m.newNode(0.0, 0.0);
	Node n1 = m.newNode(1.0, 0.0);
	Node n2 = m.newNode(1.0, 1.0);
	Node n3 = m.newNode(0.0, 1.0);
	Node n4 = m.newNode(0.5, -0.5);
	Node n5 = m.newNode(-0.5, 0.5);

	Face T1 = m.newTriangle(n0,n1,n2);
	Face T2 = m.newTriangle(n0,n2,n3);
	Face T3 = m.newTriangle(n0,n3,n5);
	Face T4 = m.newTriangle(n0,n4,n1);

	std::vector<Face> adj;
	adj.resize(3);
	adj[0] = Face();
	adj[1] = T2;
	adj[2] = T4;
	T1.set<Face>(adj);
	adj[0] =  Face();
	adj[1] = T3;
	adj[2] = T1;
	T2.set<Face>(adj);
	adj[0] =  Face();
	adj[1] =  Face();
	adj[2] = T2;
	T3.set<Face>(adj);
	adj[0] =  Face();
	adj[1] = T1;
	adj[2] = Face();
	T4.set<Face>(adj);
	bool done = algo.swap(T1,T2);

	CPPUNIT_ASSERT(done==true);
}
/*----------------------------------------------------------------------------*/
void Delaunay2DTestSuite::testSwapNonConvex(){
	std::vector<gmds::math::Point> pnts;

	Delaunay2D<gmds::TCoord,GeomToolKit >  algo(pnts);
	IGMesh& m = algo.getMesh();

	Node n0 = m.newNode(0.75, 0.75);
	Node n1 = m.newNode(1.0, 0.0);
	Node n2 = m.newNode(1.0, 1.0);
	Node n3 = m.newNode(0.0, 1.0);

	Face T1 = m.newTriangle(n0,n1,n2);
	Face T2 = m.newTriangle(n0,n2,n3);

	std::vector<Face> adj;
	adj.resize(3);
	adj[0] =  Face();
	adj[1] = T2;
	adj[2] =  Face();
	T1.set<Face>(adj);
	adj[0] =  Face();
	adj[1] =  Face();
	adj[2] = T1;
	T2.set<Face>(adj);

	bool done = algo.swap(T1,T2);
	CPPUNIT_ASSERT(done==false);
}
/*----------------------------------------------------------------------------*/
void Delaunay2DTestSuite::testGetTriangleStrip(){
	std::vector<gmds::math::Point> pnts;

	Delaunay2D<gmds::TCoord,GeomToolKit >  algo(pnts);
	IGMesh& m = algo.getMesh();

	Node n0 = m.newNode(0.0, 0.0);
	Node n1 = m.newNode(1.0, 0.0);
	Node n2 = m.newNode(1.0, 1.0);
	Node n3 = m.newNode(0.0, 1.0);

	Face T1 = m.newTriangle(n0,n1,n2);
	Face T2 = m.newTriangle(n0,n2,n3);

	std::vector<Face> adj;
	adj.resize(3);
	adj[0] = Face();
	adj[1] = T2;
	adj[2] = Face();
	T1.set<Face>(adj);
	adj[0] = Face();
	adj[1] = Face();
	adj[2] = T1;
	T2.set<Face>(adj);

	std::vector<gmds::Face> strip;
	algo.getTriangleStrip(n1,n3,strip);
	CPPUNIT_ASSERT(strip.size()==2);
	algo.getTriangleStrip(n1,n0,strip);
	CPPUNIT_ASSERT(strip.size()==0);
	algo.getTriangleStrip(n0,n2,strip);
	CPPUNIT_ASSERT(strip.size()==0);
	algo.getTriangleStrip(n2,n3,strip);
	CPPUNIT_ASSERT(strip.size()==0);
}
/*----------------------------------------------------------------------------*/
void Delaunay2DTestSuite::testBuildEdgeWith6ConvexTriangles(){
	std::vector<gmds::math::Point> pnts;

	Delaunay2D<gmds::TCoord,GeomToolKit >  algo(pnts);
	IGMesh& m = algo.getMesh();

	Node n0 = m.newNode(0.0, 0.0);
	Node n1 = m.newNode(1.0, 0.0);
	Node n2 = m.newNode(2.0, 0.0);
	Node n3 = m.newNode(3.0, 0.5);
	Node n4 = m.newNode(2.0, 1.0);
	Node n5 = m.newNode(1.0, 1.0);
	Node n6 = m.newNode(0.0, 1.0);
	Node n7 = m.newNode(-1.0, 0.5);

	Face T0 = m.newTriangle(n0,n6,n7);
	Face T1 = m.newTriangle(n6,n0,n1);
	Face T2 = m.newTriangle(n1,n2,n5);
	Face T3 = m.newTriangle(n4,n5,n2);
	Face T4 = m.newTriangle(n1,n5,n6);
	Face T5 = m.newTriangle(n4,n2,n3);

	std::vector<Face> adj;
	adj.resize(3);
	adj[0] = Face();
	adj[1] = Face();
	adj[2] = T1;
	T0.set<Face>(adj);
	adj[0] = Face();
	adj[1] = T4;
	adj[2] = T0;
	T1.set<Face>(adj);
	adj[0] = T3;
	adj[1] = T4;
	adj[2] = Face();
	T2.set<Face>(adj);
	adj[0] = T2;
	adj[1] = T5;
	adj[2] = Face();
	T3.set<Face>(adj);
	adj[0] = Face();
	adj[1] = T1;
	adj[2] = T2;
	T4.set<Face>(adj);
	adj[0] = Face();
	adj[1] = Face();
	adj[2] = T3;
	T5.set<Face>(adj);

	std::vector<gmds::Face> strip;
	algo.getTriangleStrip(n3,n7,strip);
	CPPUNIT_ASSERT(strip.size()==6);

	algo.constraintToEdge(n3,n7);
	gmds::LimaWriter<IGMesh> writer(m);
	writer.write("strip",N|F);

}
/*----------------------------------------------------------------------------*/

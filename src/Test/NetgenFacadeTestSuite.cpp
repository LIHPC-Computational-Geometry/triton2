/*----------------------------------------------------------------------------*/
/*
 * NetgenFacadeTestSuite.cpp
 *
 *  Created on: 23 mai 2011
 *      Author: ledouxf
 */
/*----------------------------------------------------------------------------*/
#include "Triton2/NetgenInterface/NetgenFacade.h"
#include "NetgenFacadeTestSuite.h"
#include <GMDS/IG/IGMesh.h>
#include <GMDSCEA/LimaWriter.h>
#include <GMDSCEA/LimaReader.h>
#include <GMDS/CAD/FacetedGeomManager.h>
#include <GMDS/Math/Point.h>
/*----------------------------------------------------------------------------*/
using namespace gmds;
using namespace triton;
using namespace nglib;
/*----------------------------------------------------------------------------*/
void NetgenFacadeTestSuite::testTet(){

	IGMesh plc(DIM3|F|N|F2N);

	Node n0 = plc.newNode(0,0,0);
	Node n1 = plc.newNode(10,0,0);
	Node n2 = plc.newNode(0,10,0);
	Node n3 = plc.newNode(0,0,10);

	IGMesh::surface& surf1 = plc.newSurface("s1");
	surf1.add(plc.newTriangle(n0,n2,n1));

	IGMesh::surface& surf2 = plc.newSurface("s2");
	surf2.add(plc.newTriangle(n0,n1,n3));

	IGMesh::surface& surf3 = plc.newSurface("s3");
	surf3.add(plc.newTriangle(n2,n0,n3));

	IGMesh::surface& surf4 = plc.newSurface("s4");
	surf4.add(plc.newTriangle(n1,n2,n3));

	NetgenFacade netgenCall;

	IGMesh volumeMesh(DIM3|R|F|F2N|N|R2N);
	netgenCall.generateTetMesh(plc,volumeMesh,1.2);

	std::cerr<<"Nb nodes: "<<volumeMesh.getNbNodes()<<std::endl;
	std::cerr<<"Nb cells: "<<volumeMesh.getNbRegions()<<std::endl;
	LimaWriter<IGMesh> w(volumeMesh);
	w.write("netgen_tet.mli",R|N);

	CPPUNIT_ASSERT(volumeMesh.getNbNodes()!=0);
	CPPUNIT_ASSERT(volumeMesh.getNbRegions()!=0);
}
/*----------------------------------------------------------------------------*/
void NetgenFacadeTestSuite::testCube(){

	IGMesh plc(DIM3|F|N|F2N);


	Node n0 = plc.newNode(0.,0.,0.);
	Node n1 = plc.newNode(0.,10.,0.);
	Node n2 = plc.newNode(10.,10.,0.);
	Node n3 = plc.newNode(10.,0.,0.);

	Node n4 = plc.newNode(0,0,10.);
	Node n5 = plc.newNode(0,10,10.);
	Node n6 = plc.newNode(10,10,10.);
	Node n7 = plc.newNode(10,0,10.);

	IGMesh::surface& surf1 = plc.newSurface("s1");
	surf1.add(plc.newTriangle(n0,n1,n2));
	surf1.add(plc.newTriangle(n0,n2,n3));

	IGMesh::surface& surf2 = plc.newSurface("s2");
	surf2.add(plc.newTriangle(n4,n6,n5));
	surf2.add(plc.newTriangle(n4,n7,n6));

	IGMesh::surface& surf3 = plc.newSurface("s3");
	surf3.add(plc.newTriangle(n0,n5,n1));
	surf3.add(plc.newTriangle(n0,n4,n5));

	IGMesh::surface& surf4 = plc.newSurface("s4");
	surf4.add(plc.newTriangle(n1,n6,n2));
	surf4.add(plc.newTriangle(n1,n5,n6));

	IGMesh::surface& surf5 = plc.newSurface("s5");
	surf5.add(plc.newTriangle(n2,n7,n3));
	surf5.add(plc.newTriangle(n2,n6,n7));

	IGMesh::surface& surf6 = plc.newSurface("s6");
	surf6.add(plc.newTriangle(n0,n7,n4));
	surf6.add(plc.newTriangle(n0,n3,n7));

	NetgenFacade netgenCall;

	IGMesh volumeMesh(DIM3|R|F|F2N|N|R2N);
	std::vector<Node> bnodes;
	IGMesh::node_iterator nit = plc.nodes_begin();
	for(;!nit.isDone();nit.next()){
		bnodes.push_back(nit.value());
	}
	std::vector<Face> bfaces;
	IGMesh::face_iterator fit = plc.faces_begin();

	for(;!fit.isDone();fit.next()){
		bfaces.push_back(fit.value());
	}
//	netgenCall.generateTetMesh(plc,volumeMesh,1.2);
	netgenCall.generateTetMesh(bnodes,bfaces,volumeMesh);


	LimaWriter<IGMesh> w(volumeMesh);
	w.write("netgen_cube.mli",R|F|N);


	CPPUNIT_ASSERT(volumeMesh.getNbNodes()!=0);
	CPPUNIT_ASSERT(volumeMesh.getNbRegions()!=0);
}
/*----------------------------------------------------------------------------*/
void NetgenFacadeTestSuite::testSurface(){

/*	geom::FacetedGeomManager m;
	m.importFAC("/cea/S/dsku/sirius/hal1/home/s3/ledouxf/workspace/GMDS/src/GMDSUnitTests/Samples/hook.fac");

	std::vector<geom::GeomSurface* > geom_surfaces;
	m.getSurfaces(geom_surfaces);

	NetgenFacade netgenCall;

	IGMesh surfaceMesh(DIM3|F|F2N|N);
	netgenCall.generateTriMesh(*(geom_surfaces[1]),m,surfaceMesh,1.2);
	VTKWriter<IGMesh> w(surfaceMesh);
	w.write("netgen_surface_facet_model.mli",F|N);

	CPPUNIT_ASSERT(surfaceMesh.getNbNodes()!=0);
	CPPUNIT_ASSERT(surfaceMesh.getNbFaces()!=0);
	*/
}
/*----------------------------------------------------------------------------*/
void NetgenFacadeTestSuite::testSTEPSolid(){


	NetgenFacade netgenCall;

	IGMesh surfaceMesh(DIM3|F|F2N|N);
	netgenCall.generateTriMeshFromSTEP(
			"/cea/S/dsku/sirius/hal1/home/s3/ledouxf/workspace/Triton/src/Test/cylindre.stp",
			surfaceMesh,
			1.2);
	LimaWriter<IGMesh> w(surfaceMesh);
	w.write("/cea/S/dsku/sirius/hal1/home/s3/ledouxf/netgen_cylindre_facet_model.mli",F|N);

	CPPUNIT_ASSERT(surfaceMesh.getNbNodes()!=0);
	CPPUNIT_ASSERT(surfaceMesh.getNbFaces()!=0);
}
/*----------------------------------------------------------------------------*/
void NetgenFacadeTestSuite::testSTEPSurface(){


	NetgenFacade netgenCall;

	IGMesh surfaceMesh(DIM3|F|F2N|N);
	netgenCall.generateTriMeshFromSTEP(
//			"/cea/S/dsku/sirius/hal1/home/s3/ledouxf/workspace/Triton/src/Test/single_surface.stp",
			"/cea/S/dsku/sirius/hal1/home/s3/ledouxf/portion_sphere.step",
//			"/cea/S/dsku/sirius/hal1/home/s3/ledouxf/quad_face.step",
			surfaceMesh,
			1.2);
	LimaWriter<IGMesh> w(surfaceMesh);
	w.write("toto.mli",F|N);

	CPPUNIT_ASSERT(surfaceMesh.getNbNodes()!=0);
	CPPUNIT_ASSERT(surfaceMesh.getNbFaces()!=0);
}
/*----------------------------------------------------------------------------*/
void NetgenFacadeTestSuite::testOCCFaceAndMeshNodes(){


	NetgenFacade netgenCall;

	IGMesh surfaceMesh(DIM3|F|F2N|N);

	// les noeuds sont ordonnes volontairement pour ce test
	std::vector<Node> corners;
	std::vector<std::vector<Node> > edges;


	std::vector<Node> e00_10, e10_11,e11_01,e01_00;
	Node n00=surfaceMesh.newNode(0, 0, 1);
	corners.push_back(n00);
	Node n10=surfaceMesh.newNode(1, 0, 1);
	corners.push_back(n10);
	Node n01=surfaceMesh.newNode(0, 1, 1);
	corners.push_back(n01);
	Node n11=surfaceMesh.newNode(1, 1, 1);
	corners.push_back(n11);

	// courbe 00 -> 10
	e00_10.push_back(n00);
	for(unsigned int i=1;i<10;i++) {
		Node n=surfaceMesh.newNode(((double)i)/10.0, 0, 1);
		e00_10.push_back(n);
	}
	e00_10.push_back(n10);

	// courbe 10 -> 11
	e10_11.push_back(n10);
	for(unsigned int i=1;i<10;i++) {
		Node n=surfaceMesh.newNode(1, ((double)i)/10.0, 1);
		i++;
		e10_11.push_back(n);
	}
	e10_11.push_back(n11);

	//courbe 11->01
	e11_01.push_back(n11);
	for(unsigned int i=1;i<10;i++) {
		Node n=surfaceMesh.newNode(1-((double)i)/10.0, 1, 1);
		e11_01.push_back(n);
	}
	e11_01.push_back(n01);
	//courbe 01 -> 00
	e01_00.push_back(n01);
	for(unsigned int i=1;i<10;i++) {
		Node n=surfaceMesh.newNode(0, 1-((double)i)/10.0, 1);
		i++;
		i++;
		e01_00.push_back(n);
	}
	e01_00.push_back(n00);

	edges.push_back(e00_10);
	edges.push_back(e10_11);
	edges.push_back(e11_01);
	edges.push_back(e01_00);

	netgenCall.generateTriMesh(
			"/cea/S/dsku/sirius/hal1/home/s3/ledouxf/quad_face.step",
			surfaceMesh, corners, edges, 1.2);
//	netgenCall.generateTriMeshFromSTEP(
//			"/cea/S/dsku/sirius/hal1/home/s3/ledouxf/quad_face.step",
//			surfaceMesh,
//			1.2);
	LimaWriter<IGMesh> w(surfaceMesh);
	w.write("toto_nodes.mli",F|N);

	CPPUNIT_ASSERT(surfaceMesh.getNbNodes()!=0);
	CPPUNIT_ASSERT(surfaceMesh.getNbFaces()!=0);
}
/*----------------------------------------------------------------------------*/
void NetgenFacadeTestSuite::testOCCSurface(){



	NetgenFacade netgenCall;

	IGMesh surfaceMesh(DIM3|F|F2N|N);

	std::vector<math::Point > pnts;
	pnts.push_back(math::Point(0.,0.,0.));
	pnts.push_back(math::Point(1.,0.,0.));
	pnts.push_back(math::Point(1.,1.,0.));
	pnts.push_back(math::Point(0.,1.,0.));

	netgenCall.generateTriMeshWithOCC(pnts,surfaceMesh,2.0);
	LimaWriter<IGMesh> w(surfaceMesh);
	w.write("/cea/S/dsku/sirius/hal1/home/s3/ledouxf/netgen_surface_facet_model.mli",F|N);

	CPPUNIT_ASSERT(surfaceMesh.getNbNodes()!=0);
	CPPUNIT_ASSERT(surfaceMesh.getNbFaces()!=0);
}
/*----------------------------------------------------------------------------*/
void NetgenFacadeTestSuite::testM3DBlock(){

/*	const int model = DIM3|R|F|E|N|R2N|F2N|E2N|R2F|F2E;
	NetgenFacade netgenCall;

	IGMesh mesh(model), m2(model);
	LimaReader<IGMesh> reader(mesh);
	reader.read("/cea/S/dsku/sirius/hal1/home/s3/ledouxf/workspace/Triton/src/Test/toto2.mli.mli",N|F);
//testPar_3
	std::vector<Node> bnodes;

	IGMesh::node_iterator itn = mesh.nodes_begin();
	for(;!itn.isDone();itn.next()){
		Node n = itn.value();

		bnodes.push_back(n);
		m2.newNode(n.X(),n.Y(),n.Z());

	}
	std::vector<Face> bfaces;
	IGMesh::face_iterator itf = mesh.faces_begin();
	for(;!itf.isDone();itf.next()){
		Face f =itf.value();
		if(f.getType()==GMDS_TRIANGLE){
			bfaces.push_back(f);

			m2.newFace(f.getIDs<Node>());

		}
	}
//	// (Z=0) NODES
//	bnodes.push_back(mesh.newNode(0.,0.,0.));
//	bnodes.push_back(mesh.newNode(1.,0.,0.));
//	bnodes.push_back(mesh.newNode(2.,0.,0.));
//	bnodes.push_back(mesh.newNode(0.,1.,0.));
//	bnodes.push_back(mesh.newNode(1.,1.,0.));
//	bnodes.push_back(mesh.newNode(2.,1.,0.));
//	bnodes.push_back(mesh.newNode(0.,2.,0.));
//	bnodes.push_back(mesh.newNode(1.,2.,0.));
//	bnodes.push_back(mesh.newNode(2.,2.,0.));
//
//	// (Z=2) NODES
//	bnodes.push_back(mesh.newNode(0.,0.,2.));
//	bnodes.push_back(mesh.newNode(1.,0.,2.));
//	bnodes.push_back(mesh.newNode(2.,0.,2.));
//	bnodes.push_back(mesh.newNode(0.,1.,2.));
//	bnodes.push_back(mesh.newNode(1.,1.,2.));
//	bnodes.push_back(mesh.newNode(2.,1.,2.));
//	bnodes.push_back(mesh.newNode(0.,2.,2.));
//	bnodes.push_back(mesh.newNode(1.,2.,2.));
//	bnodes.push_back(mesh.newNode(2.,2.,2.));
//
//	// (Z=1) NODES
//	bnodes.push_back(mesh.newNode(0.,0.,1.));
//	bnodes.push_back(mesh.newNode(1.,0.,1.));
//	bnodes.push_back(mesh.newNode(2.,0.,1.));
//	bnodes.push_back(mesh.newNode(0.,1.,1.));
//	bnodes.push_back(mesh.newNode(2.,1.,1.));
//	bnodes.push_back(mesh.newNode(0.,2.,1.));
//	bnodes.push_back(mesh.newNode(1.,2.,1.));
//	bnodes.push_back(mesh.newNode(2.,2.,1.));
//
//
//	// 1
//	bfaces.push_back(mesh.newTriangle(0,4,1));
//	bfaces.push_back(mesh.newTriangle(1,4,2));
//	bfaces.push_back(mesh.newTriangle(2,4,5));
//	bfaces.push_back(mesh.newTriangle(5,4,8));
//	bfaces.push_back(mesh.newTriangle(8,4,7));
//	bfaces.push_back(mesh.newTriangle(7,4,6));
//	bfaces.push_back(mesh.newTriangle(6,4,3));
//	bfaces.push_back(mesh.newTriangle(3,4,0));
//	// 2
//	bfaces.push_back(mesh.newTriangle(9 ,10,13));
//	bfaces.push_back(mesh.newTriangle(10,11,13));
//	bfaces.push_back(mesh.newTriangle(11,14,13));
//	bfaces.push_back(mesh.newTriangle(14,17,13));
//	bfaces.push_back(mesh.newTriangle(17,16,13));
//	bfaces.push_back(mesh.newTriangle(16,15,13));
//	bfaces.push_back(mesh.newTriangle(15,12,13));
//	bfaces.push_back(mesh.newTriangle(12,9 ,13));
//	// 3
//	bfaces.push_back(mesh.newTriangle(0 ,1 ,19));
//	bfaces.push_back(mesh.newTriangle(1 ,2 ,19));
//	bfaces.push_back(mesh.newTriangle(2 ,20,19));
//	bfaces.push_back(mesh.newTriangle(20,11,19));
//	bfaces.push_back(mesh.newTriangle(11,10,19));
//	bfaces.push_back(mesh.newTriangle(10,9 ,19));
//	bfaces.push_back(mesh.newTriangle(9 ,18,19));
//	bfaces.push_back(mesh.newTriangle(18,0 ,19));
//	// 4
//	bfaces.push_back(mesh.newTriangle(2 ,5 ,22));
//	bfaces.push_back(mesh.newTriangle(5 ,8 ,22));
//	bfaces.push_back(mesh.newTriangle(8 ,25,22));
//	bfaces.push_back(mesh.newTriangle(25,17,22));
//	bfaces.push_back(mesh.newTriangle(17,14,22));
//	bfaces.push_back(mesh.newTriangle(14,11,22));
//	bfaces.push_back(mesh.newTriangle(11,20,22));
//	bfaces.push_back(mesh.newTriangle(20,2 ,22));
//	// 5
//	bfaces.push_back(mesh.newTriangle(8 ,7 ,24));
//	bfaces.push_back(mesh.newTriangle(7 ,6 ,24));
//	bfaces.push_back(mesh.newTriangle(6 ,23,24));
//	bfaces.push_back(mesh.newTriangle(23,15,24));
//	bfaces.push_back(mesh.newTriangle(15,16,24));
//	bfaces.push_back(mesh.newTriangle(16,17,24));
//	bfaces.push_back(mesh.newTriangle(17,25,24));
//	bfaces.push_back(mesh.newTriangle(25,8 ,24));
//	// 6
//	bfaces.push_back(mesh.newTriangle(6 ,3 ,21));
//	bfaces.push_back(mesh.newTriangle(3 ,0 ,21));
//	bfaces.push_back(mesh.newTriangle(0 ,18,21));
//	bfaces.push_back(mesh.newTriangle(18,9 ,21));
//	bfaces.push_back(mesh.newTriangle(9 ,12,21));
//	bfaces.push_back(mesh.newTriangle(12,15,21));
//	bfaces.push_back(mesh.newTriangle(15,23,21));
//	bfaces.push_back(mesh.newTriangle(23,6 ,21));
	LimaWriter<IGMesh> w(m2);
	w.write("before_surf.mli",F|N);

	netgenCall.generateTetMesh(bnodes,bfaces,mesh);
	LimaWriter<IGMesh> wend(mesh);
	wend.write("after_vol.mli",R|F|N);

	CPPUNIT_ASSERT(mesh.getNbNodes()!=0);
	CPPUNIT_ASSERT(mesh.getNbFaces()!=0);*/
}
/*----------------------------------------------------------------------------*/
void NetgenFacadeTestSuite::testFacetModel(){

//
//	geom::FacetedGeomManager<TCoord> m;
//	m.importFAC("/cea/S/dsku/sirius/hal1/home/s3/ledouxf/workspace/GMDS/src/GMDSUnitTests/Samples/hook.fac");
//
//	std::vector<geom::GeomVolume<TCoord>* > geom_volumes;
//	m.getVolumes(geom_volumes);
//
//	NetgenFacade<DIM3|R|F|F2N|N|R2N> netgenCall;
//
//	Mesh<DIM3|R|F|F2N|N|R2N> volumeMesh;
//	netgenCall.generateTetMesh(*(geom_volumes[0]),m,volumeMesh,1.2);
//	VTKWriter<DIM3|R|F|F2N|N|R2N> w(volumeMesh);
//	w.write("netgen_facet_model.mli",R|F|N);
//
//	CPPUNIT_ASSERT(volumeMesh.getNbNodes()!=0);
//	CPPUNIT_ASSERT(volumeMesh.getNbRegions()!=0);
}
/*----------------------------------------------------------------------------*/

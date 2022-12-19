/*----------------------------------------------------------------------------*/
/*
 * TetgenFacadeTestSuite.cpp
 *
 *  Created on: 17 mai 2011
 *      Author: ledouxf
 */
/*----------------------------------------------------------------------------*/
#include "Triton2/TetgenInterface/TetgenFacade.h"
#include "TetgenFacadeTestSuite.h"
#include <GMDS/IG/IGMesh.h>
#include <GMDSCEA/LimaWriter.h>
#include <GMDSCEA/LimaReader.h>
/*----------------------------------------------------------------------------*/
using namespace gmds;
using namespace triton;
/*----------------------------------------------------------------------------*/
void TetgenFacadeTestSuite::testPLCCube(){

	IGMesh plc(DIM3|F|N|F2N);

	Node n0 = plc.newNode(0,0,0);
	Node n1 = plc.newNode(0,10,0);
	Node n2 = plc.newNode(10,10,0);
	Node n3 = plc.newNode(10,0,0);

	Node n4 = plc.newNode(0,0,10);
	Node n5 = plc.newNode(0,10,10);
	Node n6 = plc.newNode(10,10,10);
	Node n7 = plc.newNode(10,0,10);

	IGMesh::surface& surf1 = plc.newSurface("s1");
	surf1.add(plc.newTriangle(n0,n1,n2));
	surf1.add(plc.newTriangle(n0,n2,n3));
//	surf1.add(plc.newQuad(n0,n1,n2,n3));

	IGMesh::surface& surf2 = plc.newSurface("s2");
	surf2.add(plc.newTriangle(n4,n6,n5));
	surf2.add(plc.newTriangle(n4,n7,n6));
//	surf2.add(plc.newQuad(n4,n5,n6,n7));

	IGMesh::surface& surf3 = plc.newSurface("s3");
	surf3.add(plc.newTriangle(n0,n1,n5));
	surf3.add(plc.newTriangle(n0,n5,n4));
//	surf3.add(plc.newQuad(n0,n4,n5,n1));

	IGMesh::surface& surf4 = plc.newSurface("s4");
	surf4.add(plc.newTriangle(n1,n2,n6));
	surf4.add(plc.newTriangle(n1,n6,n5));
//	surf4.add(plc.newQuad(n1,n5,n6,n2));

	IGMesh::surface& surf5 = plc.newSurface("s5");
	surf5.add(plc.newTriangle(n2,n3,n7));
	surf5.add(plc.newTriangle(n2,n7,n6));
//	surf5.add(plc.newQuad(n2,n6,n7,n3));

	IGMesh::surface& surf6 = plc.newSurface("s6");
	surf6.add(plc.newTriangle(n3,n0,n4));
	surf6.add(plc.newTriangle(n3,n4,n7));
//	surf6.add(plc.newQuad(n3,n7,n4,n0));

	TetgenFacade tetgenCall;

	IGMesh volumeMesh(DIM3|R|F|F2N|N|R2N);
	tetgenCall.generateCDTMesh(plc,volumeMesh,1.2);

	LimaWriter<IGMesh> w(volumeMesh);
	w.write("tetgen1.mli",R|F|N);

	CPPUNIT_ASSERT(volumeMesh.getNbNodes()!=0);
	CPPUNIT_ASSERT(volumeMesh.getNbRegions()!=0);
}
/*----------------------------------------------------------------------------*/
void TetgenFacadeTestSuite::testPLCCube2(){

//	Mesh<DIM3|R|F|F2N|N|R2N> plc;
//	std::vector<Node*> bnodes;
//	Node* n0 = plc.newNode(0,0,0);
//	Node* n1 = plc.newNode(0,10,0);
//	Node* n2 = plc.newNode(10,10,0);
//	Node* n3 = plc.newNode(10,0,0);
//
//	Node* n4 = plc.newNode(0,0,10);
//	Node* n5 = plc.newNode(0,10,10);
//	Node* n6 = plc.newNode(10,10,10);
//	Node* n7 = plc.newNode(10,0,10);
//
//	bnodes.push_back(n0);
//	bnodes.push_back(n1);
//	bnodes.push_back(n2);
//	bnodes.push_back(n3);
//	bnodes.push_back(n4);
//	bnodes.push_back(n5);
//	bnodes.push_back(n6);
//	bnodes.push_back(n7);
//
//	std::vector<std::vector<Face*> > bfacets;
//
//	std::vector<Face*> f1;
//	f1.push_back(plc.newTriangle(n0,n1,n2));
//	f1.push_back(plc.newTriangle(n0,n2,n3));
//	bfacets.push_back(f1);
//
//	std::vector<Face*> f2;
//	f2.push_back(plc.newTriangle(n4,n6,n5));
//	f2.push_back(plc.newTriangle(n4,n7,n6));
//	bfacets.push_back(f2);
//
//	std::vector<Face*> f3;
//	f3.push_back(plc.newTriangle(n0,n1,n5));
//	f3.push_back(plc.newTriangle(n0,n5,n4));
//	bfacets.push_back(f3);
//
//	std::vector<Face*> f4;
//	f4.push_back(plc.newTriangle(n1,n2,n6));
//	f4.push_back(plc.newTriangle(n1,n6,n5));
//	bfacets.push_back(f4);
//
//	std::vector<Face*> f5;
//	f5.push_back(plc.newTriangle(n2,n3,n7));
//	f5.push_back(plc.newTriangle(n2,n7,n6));
//	bfacets.push_back(f5);
//
//	std::vector<Face*> f6;
//	f6.push_back(plc.newTriangle(n3,n0,n4));
//	f6.push_back(plc.newTriangle(n3,n4,n7));
//	bfacets.push_back(f6);
//
//
//	TetgenFacade<DIM3|R|F|F2N|N|R2N> tetgenCall;
//
//
//	tetgenCall.generateTetMesh(bnodes,bfacets,plc,"-Yq1.2pa0.1");
//	VTKWriter<DIM3|R|F|F2N|N|R2N> w(plc);
//	w.write("tetgen",R|F|N);
//
//	CPPUNIT_ASSERT(plc.getNbNodes()!=0);
//	CPPUNIT_ASSERT(plc.getNbRegions()!=0);
}
/*----------------------------------------------------------------------------*/
void TetgenFacadeTestSuite::testM3DBlock(){

	const int model = DIM3|R|F|E|N|R2N|F2N|E2N|R2F|F2E;

	IGMesh  mesh(model), m2(model);
	LimaReader<IGMesh> reader(mesh);
	reader.read("/cea/S/dsku/sirius/hal1/home/s3/ledouxf/workspace/Triton/src/Test/toto2.mli",N|F);
//testPar_3
	std::vector<Node> bnodes;

	IGMesh::node_iterator itn = mesh.nodes_begin();
	for(;!itn.isDone();itn.next()){
		Node n = itn.value();

		bnodes.push_back(n);
		m2.newNode(n.X(),n.Y(),n.Z());

	}

	std::vector<std::vector<Face> > bfacets;


	IGMesh::surfaces_iterator  itf = mesh.surfaces_begin();
	for(;itf!=mesh.surfaces_end();itf++){
		std::vector<Face> bfaces;
		IGMesh::surface f = *itf;
		std::vector<Face> faces = f.cells();
		for(unsigned int i_f=0;i_f<faces.size();i_f++){
			Face  ff = faces[i_f];
			if(ff.getType()==GMDS_TRIANGLE){
				bfaces.push_back(ff);

				m2.newFace(ff.getIDs<Node>());
			}
		}
		bfacets.push_back(bfaces);
	}

	LimaWriter<IGMesh> w(m2);
	w.write("TETGEN_bef_surf.mli",F|N);

	TetgenFacade tetgenCall;

	tetgenCall.generateTetMesh(bnodes,bfacets,m2,"Yq1.2pa0.1");
	LimaWriter<IGMesh> wend(m2);
	wend.write("TETGEN_after_vol.mli",R|N);

	CPPUNIT_ASSERT(mesh.getNbNodes()!=0);
	CPPUNIT_ASSERT(mesh.getNbFaces()!=0);
}
/*----------------------------------------------------------------------------*/

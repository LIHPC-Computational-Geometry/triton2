/*----------------------------------------------------------------------------*/
/*
 *  main.cpp
 *
 *  Created on: 36 mars 2013
 *      Author: ledouxf
 */
/*----------------------------------------------------------------------------*/
#include "GMDS/IG/IGMesh.h"
#include "GMDS/Math/Point.h"
#include "GMDS/IG/IGMeshDoctor.h"

#include "GMDSCEA/LimaWriter.h"
#include "GMDSCEA/LimaReader.h"
#include "GMDS/IO/VTKWriter.h"
#include "GMDS/IO/VTKReader.h"
#include "GMDS/Utils/CommonTypes.h"
//#include "GMDSMeshTools/LaplacianSmoothing.h"
//#include "GMDSMeshTools/BoundaryOperator.h"

/*----------------------------------------------------------------------------*/
#include<vector>
#include <math.h>
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
const int model = gmds::DIM2|gmds::N|gmds::E|gmds::F|gmds::F2N|gmds::N2F|gmds::E2N;
/*----------------------------------------------------------------------------*/
void addSurfaceInfo(gmds::IGMesh* m, gmds::Face& fi, gmds::Face& f){

	gmds::IGMesh::surfaces_iterator it_surf = m->surfaces_begin();
	for(;it_surf!=m->surfaces_end();it_surf++){
		gmds::IGMesh::surface& s = *it_surf;
		if(s.has(fi)){
			s.add(f);
		}
	}
}
/*----------------------------------------------------------------------------*/
void addInfo(gmds::IGMesh* m,
		gmds::Edge& e,	//arete initiale
		gmds::Edge& e1, gmds::Edge& e2, //arete decoupe de e
		gmds::Node& n12 //sommet interne ajoute
		)
{
	std::vector<gmds::Node> nodes_e = e.get<gmds::Node>();
	gmds::Node n1 = nodes_e[0];
	gmds::Node n2 = nodes_e[1];

	gmds::IGMesh::lines_iterator it_lines = m->lines_begin();
	for(;it_lines!=m->lines_end();it_lines++){
		gmds::IGMesh::line& l = *it_lines;
		if(l.has(e)){
			l.add(e1);
			l.add(e2);
		}
	}

	gmds::IGMesh::clouds_iterator it_clouds = m->clouds_begin();
	for(;it_clouds!=m->clouds_end();it_clouds++){
		gmds::IGMesh::cloud& c = *it_clouds;
		if(c.has(n1) && c.has(n2)){
			c.add(n12);
		}
	}
}
/*----------------------------------------------------------------------------*/
int main(int argc, char** argv){
	std::cout<<"==== decoupage en Quads ====\n";

	if(argc!=3){
		std::cout<<"Usage : in.mli out.mli"<<std::endl;
		return(0);
	}
	std::string file_in(argv[1]);
	std::string file_out(argv[2]);

	std::cout<<"Lecture du maillage "<<std::endl;

	gmds::IGMesh  m(model);
	gmds::LimaReader<gmds::IGMesh> reader(m);
	reader.read(file_in,gmds::F|gmds::E|gmds::N);
	gmds::IGMeshDoctor  doctor(&m);
	doctor.buildEdgesAndX2E();
	doctor.updateUpwardConnectivity();

	std::cout<<"Nb faces = "<<m.getNbFaces()<<std::endl;
	std::cout<<"Nb bras  = "<<m.getNbEdges()<<std::endl;
	std::cout<<"Nb nodes = "<<m.getNbNodes()<<std::endl;


	std::vector<gmds::Face*> toDestroy;
	std::vector<gmds::Edge*> toDestroyEdges;

	std::map<gmds::FakeEdge::ID,gmds::Node*> newNodesOnEdge;
	std::map<gmds::Face*,gmds::Node*> faceCenter;

	std::cout<<" Decoupage des bras"<<std::endl;
	gmds::Mesh<model>::edges_iterator ite = m.edges_begin();
	for(;!ite->isDone();ite->next()){
		gmds::Edge* e = ite->currentItem();
		toDestroyEdges.push_back(e);

	}

	for(unsigned int i=0;i<toDestroyEdges.size();i++){
		gmds::Edge* e = toDestroyEdges[i];

		std::vector<gmds::Node*> nodes = e->getNodes();
		Node* ni = nodes[0];
		Node* nj = nodes[1];

		double di= ni.X()*ni.X() + ni.Y()*ni.Y() + ni.Z()*ni.Z();
		double dj= nj.X()*nj.X() + nj.Y()*nj.Y() + nj.Z()*nj.Z();


		gmds::TCoord xij = (ni.X()+nj.X())/2.0;
		gmds::TCoord yij = (ni.Y()+nj.Y())/2.0;
		gmds::TCoord zij = (ni.Z()+nj.Z())/2.0;

		if(fabs(di-dj)<1e-7)
		{
			double norm_from = sqrt(xij*xij + yij *yij + zij*zij);
			double norm_to = sqrt(di);
			double x_unit = xij/norm_from;
			double y_unit = yij/norm_from;
			double z_unit = zij/norm_from;
			xij = norm_to*x_unit;
			yij = norm_to*y_unit;
			zij = norm_to*z_unit;
		}
		Node* nk = m.newNode(xij,yij,zij);

		Edge* eik = m.newEdge(ni,nk);
		Edge* ejk = m.newEdge(nj,nk);

		newNodesOnEdge[gmds::FakeEdge(ni->getID(),nj->getID()).getID()] = nk;

		addInfo(&m,e,eik,ejk,nk);
		m.deleteEdge(e);
	}

	std::cout<<" Decoupage des faces"<<std::endl;

	gmds::Mesh<model>::faces_iterator itf = m.faces_begin();

	for(;!itf->isDone();itf->next()){
		gmds::Face* f = itf->currentItem();
		toDestroy.push_back(f);

		std::vector<gmds::Node*> nodes = f->getNodes();
		int nb_nodes = nodes.size();
		//point associe a la face
		gmds::TCoord x=0, y=0, z=0;

		gmds::Node* center=0;

		bool sphericalCase=false;

		if(nb_nodes==4){
			Node n0 = nodes[0];
			Node n1 = nodes[1];
			Node n2 = nodes[2];
			Node n3 = nodes[3];

			double d0= n0.X()*n0.X() + n0.Y()*n0.Y() + n0.Z()*n0.Z();
			double d1= n1.X()*n1.X() + n1.Y()*n1.Y() + n1.Z()*n1.Z();
			double d2= n2.X()*n2.X() + n2.Y()*n2.Y() + n2.Z()*n2.Z();
			double d3= n3.X()*n3.X() + n3.Y()*n3.Y() + n3.Z()*n3.Z();

			Node nA=0, nB=0;
			if(fabs(d0-d1)<1e-7) //n0 et n1 sur le meme arc de cercle
			{
				nA = newNodesOnEdge[gmds::FakeEdge(n0->getID(),n1->getID()).getID()];
				nB = newNodesOnEdge[gmds::FakeEdge(n2->getID(),n3->getID()).getID()];
				sphericalCase=true;
			}
			else if(fabs(d0-d3)<1e-7)
			{
				nA = newNodesOnEdge[gmds::FakeEdge(n0->getID(),n3->getID()).getID()];
				nB = newNodesOnEdge[gmds::FakeEdge(n1->getID(),n2->getID()).getID()];
				sphericalCase=true;
			}
			if(sphericalCase==true)
			{
				gmds::TCoord xAB = (nA.X()+nB.X())/2.0;
				gmds::TCoord yAB = (nA.Y()+nB.Y())/2.0;
				gmds::TCoord zAB = (nA.Z()+nB.Z())/2.0;
				center = m.newNode(xAB,yAB,zAB);
			}

		}

		if(sphericalCase==false)
		{
			for(unsigned int i=0;i<nb_nodes;i++){
				Node* ni = nodes[i];
				x+=ni.X();
				y+=ni.Y();
				z+=ni.Z();
			}
			center = m.newNode(x/nb_nodes, y/nb_nodes, z/nb_nodes);
		}
		if(center==0)
			std::cout<<"ERREUR"<<std::endl;
		faceCenter[f] = center;
	}

	for(unsigned int i=0;i<toDestroy.size();i++){
		gmds::Face* fi = toDestroy[i];


		std::vector<gmds::Node*> nodes = fi->getNodes();
		int nb_nodes = nodes.size();
		gmds::Node* center = faceCenter[fi];
		//noeuds sur les aretes
		std::vector<gmds::Node*> edgeNodes;
		Node* ni = nodes[0];
		for(unsigned int j=1;j<nb_nodes;j++){
			Node* nj = nodes[j];
			gmds::FakeEdge eij(ni->getID(),nj->getID());
			gmds::Node* nij=0;
			if(newNodesOnEdge.find(eij.getID())==newNodesOnEdge.end())
			{
				//pas trouve, on le cree
				gmds::TCoord xij = (ni.X()+nj.X())/2.0;
				gmds::TCoord yij = (ni.Y()+nj.Y())/2.0;
				gmds::TCoord zij = (ni.Z()+nj.Z())/2.0;
				nij = m.newNode(xij,yij,zij);
				newNodesOnEdge[eij.getID()] = nij;
			}
			else
				nij = newNodesOnEdge[eij.getID()];

			edgeNodes.push_back(nij);
			ni=nj;
		}

		//derniere arete
		ni = nodes[nodes.size()-1];
		gmds::Node* nj = nodes[0];
		gmds::FakeEdge eij(ni->getID(),nj->getID());
		gmds::Node* nij=0;
		if(newNodesOnEdge.find(eij.getID())==newNodesOnEdge.end())
		{
			//pas trouve, on le cree
			gmds::TCoord xij = (ni.X()+nj.X())/2.0;
			gmds::TCoord yij = (ni.Y()+nj.Y())/2.0;
			gmds::TCoord zij = (ni.Z()+nj.Z())/2.0;
			nij = m.newNode(xij,yij,zij);


			newNodesOnEdge[eij.getID()] = nij;
		}
		else
			nij = newNodesOnEdge[eij.getID()];

		edgeNodes.push_back(nij);

		//creation des faces
		for(unsigned int j=0;j<nb_nodes;j++){
			Node* nj = nodes[j];
			Node* ni = (j!=0)?nodes[j-1]:nodes[nodes.size()-1];
			Node* nk = (j!=nodes.size()-1)?nodes[j+1]:nodes[0];
			gmds::FakeEdge eij(ni->getID(),nj->getID());
			gmds::FakeEdge ejk(nj->getID(),nk->getID());
			gmds::Node* nij = newNodesOnEdge[eij.getID()];
			gmds::Node* njk = newNodesOnEdge[ejk.getID()];
			gmds::Face *q = m.newQuad(nij,nj,njk,center);

			addSurfaceInfo(&m,fi,q);

		}

		m.deleteFace(fi);

	}
	std::cout<<"Ecriture du maillage"<<std::endl;
	gmds::LimaWriter<model> writer(m);
	writer.write(file_out,N|E|F);
	std::cout<<"--> done\n";

	return 0;
}
/*----------------------------------------------------------------------------*/

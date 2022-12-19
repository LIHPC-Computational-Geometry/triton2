/*----------------------------------------------------------------------------*/
/*
 * TritonPluginAPI.cpp
 *
 *  Created on: 24 avril 2013
 *      Author: ledouxf
 */
/*----------------------------------------------------------------------------*/
#include "TritonPluginAPI.h"
#include "Triton2/Core/Delaunay2D.h"
/*----------------------------------------------------------------------------*/
#include <string>
#include <iostream>
/*----------------------------------------------------------------------------*/
using namespace triton;
/*----------------------------------------------------------------------------*/
TritonPluginAPI::TritonPluginAPI(){;}
/*----------------------------------------------------------------------------*/
TritonPluginAPI::~TritonPluginAPI(){;}
/*----------------------------------------------------------------------------*/
void TritonPluginAPI::triangulate(
		const std::vector<double>& ANodeX,
		const std::vector<double>& ANodeY,
		/** (I/O) Nombre de mailles */
		size_t&         nb_mesh,
		/** (I/O) Nombre total de noeuds (interieur + contour) */
		size_t&         nb_node,
		/** (I/O) Nombre de noeuds internes */
		size_t&         nb_int_node,
		/** (I/O) Nombre de noeuds du contour */
		size_t&         nb_out_node,
		/** (I/O) Coordonnees X des noeuds internes */
		double        *& int_node_x,
		/** (I/O) Coordonnees Y des noeuds internes */
		double        *& int_node_y,
		/** (I/O) Tableau du nombre de noeuds par maille */
		size_t        *& nb_node_per_mesh,
		/** (I/O) Table de connectivite noeuds/mailles */
		size_t        *& mesh_to_node)
{
	std::cout<<"Triangulation avec Triton"<<std::endl;
	std::vector<math::Point > boundaryNodes;


	for(unsigned int i=0;i<ANodeX.size();i++){
		boundaryNodes.push_back(math::Point(ANodeX[i], ANodeY[i]) );
	}

	Delaunay2D<gmds::TCoord,GeomToolKit >::OrientedBoundary bnd(boundaryNodes);

	Delaunay2D<gmds::TCoord,GeomToolKit > algo(bnd);
	algo.mesh();
	gmds::Mesh<DIM2|N|F|F2N|F2F>& m = algo.getMesh();

	std::cout<<"Nb mailles = "<<m.getNbFaces()<<std::endl;
	std::cout<<"Nb noeuds  = "<<m.getNbNodes()<<std::endl;

	nb_node = m.getNbNodes();
	nb_out_node = boundaryNodes.size();
	nb_int_node = nb_node - nb_out_node;
	std::cout<<"(T, O, I) = ("<<nb_node<<", "<<nb_out_node<<", "<<nb_int_node<<")"<<std::endl;
	int_node_x = new double[nb_int_node];
	int_node_y = new double[nb_int_node];

	std::map<gmds::id,int> gmdsIDs_toMgxID;

	gmds::Mesh<DIM2|N|F|F2N|F2F>::nodes_iterator it_node     = m.nodes_begin();
	// la numerotation commence
	int i_in=0;
	int i_total=0;
	int bmark = algo.getBoundaryMark();
	for(;!it_node->isDone();it_node->next()){
		Node *n = it_node->currentItem();
		if(!m.isMarked(n,bmark)){
			int_node_x[i_in] = n->getX().toDouble();
			int_node_y[i_in] = n->getY().toDouble();
			gmdsIDs_toMgxID[n->getID()]=nb_out_node+i_in;
			i_in++;
			i_total++;
		}
		else{

			gmdsIDs_toMgxID[n->getID()]= i_total;
			i_total++;
		}
	}

	nb_mesh=m.getNbFaces();
	nb_node_per_mesh= new size_t[nb_mesh];
	for(unsigned int i=0;i<nb_mesh;i++)
		nb_node_per_mesh[i]=3;

	gmds::Mesh<DIM2|N|F|F2N|F2F>::faces_iterator it_face = m.faces_begin();
	mesh_to_node=new size_t[3*nb_mesh];
	int j=0;
	for(;!it_face->isDone();it_face->next()){
		Face *f = it_face->currentItem();
		std::vector<gmds::id> nodeIDs = f->getNodeIDs();
		mesh_to_node[3*j]   = gmdsIDs_toMgxID[nodeIDs[0]];
		mesh_to_node[3*j+1] = gmdsIDs_toMgxID[nodeIDs[1]];
		mesh_to_node[3*j+2] = gmdsIDs_toMgxID[nodeIDs[2]];
		j++;
	}

}
/*----------------------------------------------------------------------------*/

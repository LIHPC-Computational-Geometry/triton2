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
#include "GMDSCEA/LimaWriter.h"
#include "GMDSCEA/LimaReader.h"
#include "GMDS/IO/VTKWriter.h"
#include "GMDS/IO/VTKReader.h"
#include "GMDSMeshTools/LaplacianSmoothing.h"
#include "GMDSMeshTools/BoundaryOperator.h"
#include "GMDSMeshTools/SheetOperator.h"
#include "Triton2/Core/Delaunay2D.h"
/*----------------------------------------------------------------------------*/
#include<vector>
#include <sstream>
#include <math.h>
/*----------------------------------------------------------------------------*/
using namespace triton;
/*----------------------------------------------------------------------------*/
const int mesh_model = DIM2|N|F|F2N|N2F;
/*----------------------------------------------------------------------------*/
std::set<Node*> getAdjNodes(Node* n){
	std::set<Node*> adj_nodes;
	std::vector<Face*> n_faces = n->getFaces();
	for(unsigned int i_face=0;i_face<n_faces.size();i_face++){
		Face *current_face = n_faces[i_face];
		Node *n1, *n2;
		current_face->getAdjacentNodes(n,&n1,&n2);
		adj_nodes.insert(n1);
		adj_nodes.insert(n2);
	}

	return adj_nodes;
}
/*----------------------------------------------------------------------------*/
void randomInit()
{

  unsigned int seed;

  std::string input_file_name;
  input_file_name.append("/dev/random");
  std::ifstream input_file_stream(input_file_name.c_str(),std::ios::in);
  input_file_stream >> seed;
  input_file_stream.close();

  //seed = 2;
  srand(time(NULL) + seed);
}
/*----------------------------------------------------------------------------*/
void computeVolumes(gmds::Mesh<mesh_model> &m,
		double& v_leger, double&v_lourd, int mark_leger, int mark_lourd)
{
	double pi_div_3 = 3.1415926535897931/3.0;
	gmds::Mesh<mesh_model>::faces_iterator it_f= m.faces_begin();
	for(;!it_f->isDone();it_f->next())
	{
		Face *f = (*it_f).currentItem();
		std::vector<Node*> nodes = f->getNodes();
		int nb_nodes = nodes.size();
		double x1,y1,x2,y2;
		for(unsigned int i=0;i<nb_nodes;i++){
			Node *n1 = nodes[i];
			Node *n2 = nodes[(i+1)%nb_nodes];
    		x1 = n1->getX().toDouble();
    		x2 = n2->getX().toDouble();
    		y1 = n1->getY().toDouble();
    		y2 = n2->getY().toDouble();

		}
		double vol =(x2-x1)*(y2*y2+y1*y2+y1*y1);

		vol*= pi_div_3;
		if(vol<0)
			vol = -vol;
		if(m.isMarked(f,mark_lourd))
			v_lourd += vol;
		else if(m.isMarked(f,mark_leger))
			v_leger += vol;

	}

}
/*----------------------------------------------------------------------------*/
bool selectPath(gmds::Mesh<mesh_model> &m, gmds::Node* n, int boundary, std::vector<Node*>& nodes)
{

	//we get all the adjacent nodes of n that are inside the mesh
	std::set<Node*> adj_nodes = getAdjNodes(n);
	int nb_adj=0;
	std::set<Node*>::iterator it_nodes;

	Node *n0 = 0,  *n1 = 0;
	bool found = false;
	for(it_nodes = adj_nodes.begin(); it_nodes!=adj_nodes.end() && !found;it_nodes++){
		Node *current = *it_nodes;
		if(!m.isMarked(current,boundary)){
			n0 = current;
			found = true;
		}
	}
	found =false;
	for(it_nodes = adj_nodes.begin(); it_nodes!=adj_nodes.end() && !found;it_nodes++){
		Node *current = *it_nodes;
		if(!m.isMarked(current,boundary) && current!=n0){
			std::vector<id> fn0 = n0->getFaceIDs();
			std::vector<id> fn1 = current->getFaceIDs();
			bool sameFace = false;
			for(unsigned int face_index=0; face_index<fn0.size() && !sameFace;face_index++){
				if(std::find(fn1.begin(),fn1.end(),fn0[face_index])!=fn1.end())
					sameFace = true;
			}
			if(!sameFace){
				n1 = current;
				found = true;
			}
		}
	}

	//we have our two adjacent nodes

	if(n0==0 || n1==0){
		return false ;
	}
	Node *n0prev=0, *n1next=0;

	adj_nodes = getAdjNodes(n0);
	bool found_next = false;
	for(it_nodes = adj_nodes.begin(); it_nodes!=adj_nodes.end() && !found_next;it_nodes++)
	{
		Node *current = *it_nodes;
		if(current!=n && current!=n1 && current!=n0 && !m.isMarked(current,boundary)){
			n0prev = current;
			found_next=true;
		}
	}
	adj_nodes = getAdjNodes(n1);

	found_next = false;
	for(it_nodes = adj_nodes.begin(); it_nodes!=adj_nodes.end() && !found_next;it_nodes++)
	{
		Node *current = *it_nodes;
		if(current!=n && current!=n1 && current!=n0 && current!=n0prev  && !m.isMarked(current,boundary)){
			n1next = current;
			found_next=true;
		}
	}

	if(n0==0 || n1==0 || n0prev==0 || n1next==0){
		return false;
	}

	//sheet definition
	nodes.clear();
	nodes.push_back(n0prev);
	nodes.push_back(n0);
	nodes.push_back(n);
	nodes.push_back(n1);
	nodes.push_back(n1next);

	return true;
}
/*----------------------------------------------------------------------------*/
int main(int argc, char** argv){
	std::cout<<"==== TRI QUAD LAYERS ====\n";

	if(argc!=6){
		std::cout<<"Usage : insertLayer in.mli mat densite_leger densite_lourd densite_moy "<<std::endl;
		std::cout<<"- in.mli \t fichier de maillage a modifier"<<std::endl;
		std::cout<<"- mat\t \t nom de la zone dans laquelle ajouter des mailles fines"<<std::endl;
		std::cout<<"- densite_leger\t densite des mailles de la zone \"mat\""<<std::endl;
		std::cout<<"- densite_lourd\t densite des mailles a ajouter "<<std::endl;
		std::cout<<"- densite_moy\t densite moyenne a obtenir dans la zone \"mat\" apres insertion"<<std::endl;
		return(0);
	}


	std::string in_file = (argv[1]);
	std::string mat_name = (argv[2]);
	double densite_leger = atof(argv[3]);
	double densite_lourd = atof(argv[4]);
	double densite_moyen = atof(argv[5]);


	std::cout<<"1 - Lecture du maillage "<<std::endl;
	gmds::Mesh<mesh_model> m;
	gmds::LimaReader<mesh_model> reader(m);
	reader.read(in_file,N|F);

	gmds::MeshDoctor<mesh_model> doctor2(m);
	doctor2.updateUpwardConnectivity();



	std::cout<<"2 - recuperation de la surface "<<std::endl;
	gmds::Mesh<mesh_model>::surface &surf_mat = m.getSurface(mat_name);
	int markLight = m.getNewMark();
	std::vector<Face*>& mat_faces =  surf_mat.cells();
	for(unsigned int i=0;i<mat_faces.size();i++){
		m.mark(mat_faces[i],markLight);
	}

	int boundary = m.getNewMark();
	gmds::BoundaryOperator<mesh_model> bound_op(m);
	bound_op.markBoundaryNodes(boundary);
	gmds::Mesh<mesh_model>::faces_iterator it_f= m.faces_begin();
	for(;!it_f->isDone();it_f->next())
	{
		Face *f = (*it_f).currentItem();
		if(!m.isMarked(f,markLight)){
			std::vector<Node*> fnodes = f->getNodes();
			for(unsigned int i=0;i<fnodes.size();i++)
				m.mark(fnodes[i],boundary);
		}
	}

	std::cout<<"3 - insertion des couches fines"<<std::endl;


	//we mark all the boundary nodes to avoid to select them for
	// the sheet insertion process
	SheetOperator<mesh_model> sheet_inserter(m);

	gmds::Mesh<mesh_model>::surface &s = m.newSurface("ADD");

	randomInit();
	bool doInsertion = true;
	int nb_insertions = 0;

	//marque pour indiquer quelles faces appartiennent au leger
	int markHeavy = m.getNewMark();
	while (doInsertion){

		int nbNodes = m.getNbNodes();
		int random_val = rand();
		id id_node =(int)(((double)random_val/(double)RAND_MAX)*(double)nbNodes);

		if(id_node==0)
			id_node= 12;
		Node* n = m.getLNode(id_node);

		//je ne traite pas les sommets du bord
		if(m.isMarked(n,boundary)){
			continue;
		}
		//ni ceux qui ne sont pas dans la surface
		std::vector<Face*> n_faces = n->getFaces();
		bool is_in=false;
		for(unsigned int iface=0; iface<n_faces.size() & !is_in; iface++)
		{
			if(surf_mat.has(n_faces[iface]))
				is_in=true;
		}
		if(!is_in)
			continue;

		//======================================================
		// SELECTION DES NOEUDS DU CHEMIN A OUVRIR
		//======================================================
		std::vector<Node*> nodes;
		bool findPath = selectPath(m,n,boundary,nodes);

		if(!findPath)
			continue;
		//======================================================
		// AJOUT DE FACES QUAD ET TRIANGULAIRES
		//======================================================
		std::vector<Face*> faces = sheet_inserter.inflate(nodes);
//
//		for(unsigned int iNodes=0;iNodes<nodes.size();iNodes++)
//			m2.mark(nodes[iNodes],boundary);

		//marquage de toutes les nouvelles faces
		int markSheet = m.getNewMark();
		for(unsigned int iFace=0;iFace<faces.size();iFace++){
					m.mark(faces[iFace],markSheet);
					m.mark(faces[iFace],markHeavy);
					s.add(faces[iFace]);
					surf_mat.add(faces[iFace]);
		}

		std::set<Node*> sheetNodes;
		for(unsigned int iFace=0;iFace<faces.size();iFace++){
			std::vector<Node*> currentNodes = faces[iFace]->getNodes();
			sheetNodes.insert(currentNodes.begin(),currentNodes.end());

			//pour ne pas reutiliser de nouveaux noeuds
			for(unsigned int iNode =0;iNode<currentNodes.size();iNode++)
				m.mark(currentNodes[iNode],boundary);

		}
		std::set<Node*>::iterator it_sheetNodes=sheetNodes.begin();
		while(it_sheetNodes!=sheetNodes.end()){
			m.mark(*it_sheetNodes,markSheet);
			it_sheetNodes++;
		}

		/* Pour chaque noeud, on recupere les noeuds adjacents qui ne sont pas eux
		 * meme au bord de faces */
		it_sheetNodes=sheetNodes.begin();
		while(it_sheetNodes!=sheetNodes.end()){
			Node* currentN = *it_sheetNodes;
			std::vector<Face*> currentFaces = currentN->getFaces();
			int nbSheetFaces=0;
			for(unsigned int iF=0;iF<currentFaces.size();iF++){
				if(m.isMarked(currentFaces[iF],markSheet))
					nbSheetFaces++;
			}
			if(nbSheetFaces==2){ //noeud pas à l'extrémité de la chaine
				//on reparcourt les voisins pour déterminer un déplacement
				std::vector<TCoord> xs;
				std::vector<TCoord> ys;
				TCoord nx=currentN->getX();
				TCoord ny=currentN->getY();
				for(unsigned int iF=0;iF<currentFaces.size();iF++){
					if(!m.isMarked(currentFaces[iF],markSheet))
					{
						std::vector<Node*> nodesF=currentFaces[iF]->getNodes();
						for(unsigned int iN=0; iN<nodesF.size();iN++){
							Node *ni=nodesF[iN];
							if(!m.isMarked(ni,markSheet)){
								xs.push_back(ni->getX()-nx);
								ys.push_back(ni->getY()-ny);
							}
						}
					}
				}

				TCoord vx=0, vy=0;
				for(unsigned int iVect=0; iVect<xs.size(); iVect++){
					vx+=xs[iVect];
					vy+=ys[iVect];
				}
				if(!xs.empty()){
					vx = vx/xs.size();
					vy = vy/ys.size();
//					std::cout<<nx<<" "<<ny<<" -- "<<vx<<" "<<vy<<std::endl;
					currentN->setX(nx+(vx*0.1));
					currentN->setY(ny+(vy*0.1));
				}
			}
			//noeud suivant
			it_sheetNodes++;
		}

		//demarquage de toutes les nouvelles faces
		for(unsigned int iFace=0;iFace<faces.size();iFace++){
					m.unmark(faces[iFace],markSheet);
		}
		it_sheetNodes=sheetNodes.begin();
		while(it_sheetNodes!=sheetNodes.end()){
			m.unmark(*it_sheetNodes, markSheet);
			it_sheetNodes++;
		}
		m.freeMark(markSheet);

		double vol_leger=0;
		double vol_lourd=0;
		computeVolumes(m,vol_leger,vol_lourd, markLight,markHeavy);


		double vol_total = vol_leger+vol_lourd;

		double densite_moyen_tmp = (densite_leger*vol_leger + densite_lourd*vol_lourd)/vol_total;
//
//			gmds::LaplacianSmoothing<mesh_model> smooth_op(m);
//			smooth_op.perform(1,boundary);
//			std::cout<<"smoothing done"<<std::endl;


//		std::stringstream filename;
//		filename<<"mixed_mesh_"<<nb_insertions;
//		gmds::VTKWriter<DIM2|N|F|F2N|N2F> writer2(m);
//		writer2.write(filename.str(),N|F);
		nb_insertions++;

		double pourcentage = densite_moyen_tmp/densite_moyen;

		if(nb_insertions%50==0){
			std::cout<<"   > couches fines insérées : "<<nb_insertions;
			std::cout<<"  -- densité moyenne actuelle : "<<densite_moyen_tmp<<std::endl;
			std::cout<<nb_insertions<<" - (vol leger, vol lourd) = ("<<vol_leger<<", "<<vol_lourd<<")";
		}
		if(pourcentage> .95){
			doInsertion=false;
			std::cout<<nb_insertions<<" - (vol leger, vol lourd) = ("<<vol_leger<<", "<<vol_lourd<<")";
			std::cout<<" -> "<<densite_moyen<<" <> "<<densite_moyen_tmp<<std::endl;
			std::cout<<" -> "<<pourcentage*100<<std::endl;
		}
	}

	m.unmarkAll(boundary);
	m.freeMark(boundary);

		//mesh topology correction
//		doctor2.buildEdgesAndX2E();
//		doctor2.updateUpwardConnectivity();

	std::cout<<"STEP 3 -  MESH WRITING "<<std::endl;
	gmds::VTKWriter<DIM2|N|F|F2N|N2F> writer2(m);
	writer2.write("mixed_mesh",N|F);

	gmds::LimaWriter<DIM2|N|F|F2N|N2F> writerLima(m);
	writerLima.write("mixed_mesh.unf",N|F);

	std::cout<<"VTK Writing done\n";
	return 0;
}
/*----------------------------------------------------------------------------*/

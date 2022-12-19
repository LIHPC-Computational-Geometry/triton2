
#include "TritonPluginAPI.h"
#include <iostream>
#include <vector>
/*----------------------------------------------------------------------------*/
#include "MgxAlgo.h"

int insertionDeFailles(T_MagixMeshAlgoData *data)
{
	return 0;
}

int maillageTriangulaire(T_MagixMeshAlgoData *data)
{
	int error_code;
	int nbIter = MgxMeshAlgoGetInteger(data,"niter",&error_code);

	/** (I) Nombre total de noeuds du contour */
	int nbContours = data->nb_outline;
	std::cout<<"nb de contours: "<<nbContours<<std::endl;

	size_t* nbNodesPerSide = data->nb_node_per_outline;

	for(unsigned int i=0; i<nbContours;i++)
		std::cout<<"nb noeuds sur contour "<<i<<": "<<nbNodesPerSide[i]<<std::endl;

	std::cout<<"Maillage (0) ou lissage (1): "<<data->is_smooth<<std::endl;

	double *nodesX=data->out_node_x;
	double *nodesY=data->out_node_y;

	double *ptrX = &nodesX[0];
	double *ptrY = &nodesY[0];

	std::vector<double> x_nodes, y_nodes;

	for(unsigned int i=0; i<nbContours;i++){
		int nbNodes = nbNodesPerSide[i];
		std::cout<<"Contour "<<i<<": "<<nbNodes<<std::endl;
		for(unsigned int j=0;j<nbNodes-1;j++){
			std::cout<<"("<<*ptrX<<", "<<*ptrY<<") ";
			x_nodes.push_back(*ptrX);
			y_nodes.push_back(*ptrY);
			if(j<nbNodes-1){
				ptrX++; ptrY++;
			}
		}
		std::cout<<std::endl;
	}

	for(unsigned int i=0;i<x_nodes.size();i++){
		std::cout<<x_nodes[i]<<", "<<y_nodes[i]<<std::endl;
	}
	std::cout<<"--> nb nodes = "<<x_nodes.size()<<std::endl;

	triton::TritonPluginAPI plugin_object;
	plugin_object.triangulate(x_nodes, y_nodes,
			data->nb_mesh,
			data->nb_node,
			data->nb_int_node,
			data->nb_out_node,
			data->int_node_x,
			data->int_node_y,
			data->nb_node_per_mesh,
			data->mesh_to_node);

	  /** (I/O) Nombre de mailles */
	  size_t         nb_mesh;
	  /** (I/O) Nombre total de noeuds (interieur + contour) */
	  size_t         nb_node;
	  /** (I/O) Nombre de noeuds internes */
	  size_t         nb_int_node;
	  /** (I/O) Nombre de noeuds du contour */
	  size_t         nb_out_node;
	  /** (I/O) Coordonnees X des noeuds internes */
	  double        *int_node_x;
	  /** (I/O) Coordonnees Y des noeuds internes */
	  double        *int_node_y;
	  /** (I/O) Tableau du nombre de noeuds par maille */
	  size_t        *nb_node_per_mesh;
	  /** (I/O) Table de connectivite noeuds/mailles */
	  size_t        *mesh_to_node;
	return 0;


}

/*----------------------------------------------------------------------------*/
/*
 * TritonPluginAPI.h
 *
 *  Created on: 24 avril 2013
 *      Author: ledouxf
 */
/*----------------------------------------------------------------------------*/
#ifndef TRITON_PLUGIN_API_H_
#define TRITON_PLUGIN_API_H_
/*----------------------------------------------------------------------------*/
#include "MgxAlgo.h"
#include <vector>
/*----------------------------------------------------------------------------*/
extern "C"{

	int insertionDeFailles(T_MagixMeshAlgoData *data);

	int maillageTriangulaire(T_MagixMeshAlgoData *data);
}
namespace triton{
/*----------------------------------------------------------------------------*/
/** \class TritonPluginAPI
 *  \brief ...
 */
class TritonPluginAPI{

public:

	/*------------------------------------------------------------------------*/
    /** \brief  Constructor.
     *
     *	\param ..
     */
	TritonPluginAPI();

	/*------------------------------------------------------------------------*/
    /** \brief  Destructor.
	*/
	virtual ~TritonPluginAPI();

	/*------------------------------------------------------------------------*/
	/** \brief  Generation of a Triangular mesh from boundary nodes.
	 */
	void triangulate(const std::vector<double>& ANodeX,
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
			size_t        *& mesh_to_node);
};
/*----------------------------------------------------------------------------*/
}
/*----------------------------------------------------------------------------*/
#endif /* TRITON_PLUGIN_API_H_ */
/*----------------------------------------------------------------------------*/

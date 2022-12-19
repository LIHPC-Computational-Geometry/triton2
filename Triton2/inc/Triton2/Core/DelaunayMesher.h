/*----------------------------------------------------------------------------*/
/*
 * DelaunayMesher.h
 *
 *  Created on: 22 mai 2013
 *      Author: ledouxf
 */
/*----------------------------------------------------------------------------*/
#ifndef DELAUNAYMESHER_H_
#define DELAUNAYMESHER_H_
/*----------------------------------------------------------------------------*/
#include "GMDSMesh/Mesh.h"
#include "Gepeto/Point.h"
#include "GMDSMesh/MeshDoctor.h"
#include "GeomToolKit.h"
#include "GMDSCommon/Timer.h"
#include "GMDSIO/VTKWriter.h"
#include <iostream>
#include <sstream>
#include <string>
#include <math.h>
#include <map>
#include <list>
/*----------------------------------------------------------------------------*/
#include <vector>
/*----------------------------------------------------------------------------*/
namespace triton{
/*----------------------------------------------------------------------------*/
/** \class DelaunayMesher
 *  \brief Implements a Delaunay Algorithm starting from a set of nodes + a
 *  	   set of edges.
 *  \param TGeomKernel the class TGeomKernel must provide a geometric kernel
 *  	   responsible of computing some geometric predicates/operations. In
 *  	   order to be compliant with gmds, it must use TNum as
 *  	   numeric basic type.
 */
template <	typename TNum,
			template<typename T=TNum> class TGeomKernel>
class DelaunayMesher{

public:
	/*------------------------------------------------------------------------*/
    /** \brief  Constructor.
     *
     *	\param APoints Points we want to build a triangulation on.
     */
	DelaunayMesher(const std::vector<GEPETO::Point<2,TNum> >& APoints);

	/*------------------------------------------------------------------------*/
    /** \brief  Destructor.
	*/
	virtual ~DelaunayMesher();

	/*------------------------------------------------------------------------*/
    /** \brief  Build a triangulation on the points stored in (*this)
	*/
	void triangulate();

	/*------------------------------------------------------------------------*/
    /** \brief  Build a constrained triangulation on the points and segments
     * 			stored in (*this)
	*/
	void triangulateConstrained();

	/*------------------------------------------------------------------------*/
    /** \brief  Build a mesh on the points and segments stored in (*this)
	*/
	void mesh();
	/*------------------------------------------------------------------------*/
    /** \brief  return a reference on the resulting mesh
	*/
	gmds::Mesh<DIM2|N|F|F2N|F2F|N2F>& getMesh() {return mesh_;}

	/*------------------------------------------------------------------------*/
    /** \brief  swap faces AF1 and AF2 along their common edge. If they don't
     * 			share an edge, the swap does nothing. The creation of inverted
     * 			faces is allowed if AB is true, otherwise the swap is not
     * 			applied.
     *
     * 			To run this operation, it is assumed that faces AF1 and AF2 are
     * 			triangles and oriented in the same way.
     *
     * 	/return true if the swap was done, false otherwise
     *
	 */
	bool swap(gmds::Face*& AF1, gmds::Face*& AF2, const bool& AB = false);

	/*------------------------------------------------------------------------*/
    /** \brief  provide the strip of faces intersected by the segment [AN1,AN2]
     *
     *
	 */
	void getTriangleStrip(	gmds::Node*& AN1, gmds::Node*& AN2,
							std::vector<gmds::Face*>& AStrip);

	/*------------------------------------------------------------------------*/
    /** \brief  constraint the triangular mesh to own the edge forming by nodes
     * 			AN1 and AN2
     *
     *
	 */
	void constraintToEdge(gmds::Node*& AN1, gmds::Node*& AN2);

	int getBoundaryMark() const {return init_mark_;}

protected:

	void middleEdgeInsertion();

	/*------------------------------------------------------------------------*/
    /** \brief  class used to order 2D points in a std::set
	*/
	class PntComp{

	public:

		int operator()(const GEPETO::Point<2,TNum>& p1,
		const GEPETO::Point<2,TNum>& p2) const {
			return ((p1.getX()<p2.getX()) || ((p1.getX()==p2.getX()) && (p1.getY()<p2.getY())) );}
	};

	/*------------------------------------------------------------------------*/
    /** \brief  Build a triangulation on the points stored in (*this)
	*/
	void initTriangulation(GEPETO::Point<2,TNum>& AP_xmin_ymin,
			GEPETO::Point<2,TNum>& AP_xmax_ymin,
			GEPETO::Point<2,TNum>& AP_xmax_ymax,
			GEPETO::Point<2,TNum>& AP_xmin_ymax);

	bool isInCircumCircle(gmds::Face* AFace,
						  GEPETO::Point<2,TNum>& p);

	void buildCavity(gmds::Face* AFace, GEPETO::Point<2,TNum>& p,
					 std::vector<gmds::Face*>& cavity,
					 const bool& relax=false);

	void splitCavity(gmds::Node* ANodeToInsert,
					 std::vector<gmds::Face*>& cavity);

	gmds::Face* findFirstTriangle(GEPETO::Point<2,TNum>& p);

	void insertNode(gmds::Node* ANode);

	void removeOuterTriangles();
	/*------------------------------------------------------------------------*/
    /** \brief  Starting from an existing triangulation, this operation
     * 			constraint the triangulation to capture the edges of
     * 			init_segments_.
	*/
	void constraintBoundary();

	/*------------------------------------------------------------------------*/
    /** \brief  Check if the edge [AN1,AN2] is available in the triangulation
	*/
	bool checkEdgeExistence(Node* AN1, Node* AN2);

	/*------------------------------------------------------------------------*/
    /** \brief  Check if the point (x,y) is out of the domain
	*/
	bool isOutOfTheDomain(const TNum& AX, const TNum& AY);

	void displayMeshInfo();

	bool isBoundaryEdge(Node* AN1, Node* AN2);

	void computeCircumCenter(Face* f, GEPETO::Point<2,TNum>& center);

private:

	int init_mark_;

	std::vector<std::pair<Node*, Node* > > init_segments_;
	std::map<id,std::vector<id> > init_boundary_association_;

	std::vector<Node*> bounding_box_nodes_;
	std::vector<Node*> init_nodes_;
	std::vector<GEPETO::Point<2,TNum> > init_points_;


	gmds::Mesh<DIM2|N|F|F2N|F2F|N2F> mesh_;
};
/*----------------------------------------------------------------------------*/
#include "DelaunayMesher.t.h"
/*----------------------------------------------------------------------------*/
}
/*----------------------------------------------------------------------------*/
#endif /* DELAUNAYMESHER_H_ */
/*----------------------------------------------------------------------------*/

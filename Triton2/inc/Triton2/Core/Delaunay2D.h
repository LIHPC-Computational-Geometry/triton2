/*----------------------------------------------------------------------------*/
/*
 * Delaunay2D.h
 *
 *  Created on: 19 janv. 2011
 *      Author: ledouxf
 */
/*----------------------------------------------------------------------------*/
#ifndef DELAUNAY2D_H_
#define DELAUNAY2D_H_
/*----------------------------------------------------------------------------*/
#include "GMDS/IG/IGMesh.h"
#include "GMDS/Math/Point.h"
#include "GMDS/IG/IGMeshDoctor.h"
#include "GMDSCEA/LimaWriter.h"
#include "GeomToolKit.h"
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
/** \class Delaunay2D
 *  \brief Implements a 2D Delaunay Algorithm starting from a set of nodes + a
 *  	   set of edges.
 *  \param TGeomKernel the class TGeomKernel must provide a geometric kernel
 *  	   responsible of computing some geometric predicates/operations. In
 *  	   order to be compliant with gmds, it must use TNum as
 *  	   numeric basic type.
 */
template <	typename TNum,
			template<typename T=TNum> class TGeomKernel>
class Delaunay2D{

public:

	enum TAlgoStrategy{
		DELAUNAY_STRATEGY,
		EDGE_LENGTH_STRATEGY
	} ;
	/*------------------------------------------------------------------------*/
    /** \class OrientedBoundary
     *	\brief defines how to get a boundary to be meshed. In other word, when
     *		   you want to mesh a 2D surface whose the boundary is already
     *		   discretized, you must provide to the algorithm the ordered list
     *		   of the boundary points by filling an instance of the
     *		   OrientedBoundary class.
     *
     *		   We assumet there exists an edge linking the first point and the
     *		   last point in APoints.
     */
	class OrientedBoundary{
	public:
		OrientedBoundary(const std::vector<gmds::math::Point >& APoints);
		const std::vector<gmds::math::Point >& getOrderedPoints() const ;

	private:
		const std::vector<gmds::math::Point >& points_;
	};

	/*------------------------------------------------------------------------*/
    /** \brief  Constructor.
     *
     *	\param APoints Points we want to build a triangulation on.
     */
	Delaunay2D(const std::vector<gmds::math::Point >& APoints);

	/*------------------------------------------------------------------------*/
    /** \brief  Constructor.
     *
     *	\param ABoundary the outter boundary of the domain we want to mesh.
     */
	Delaunay2D(const OrientedBoundary& ABoundary);

	/*------------------------------------------------------------------------*/
    /** \brief  Constructor.
     *
     *	\param ABoundary the outter boundary of the domain we want to mesh.
     */
	Delaunay2D(const OrientedBoundary& AOutBoundary,
			const OrientedBoundary* AInBoundary, const int ANbIn);

	/*------------------------------------------------------------------------*/
    /** \brief  Destructor.
	*/
	virtual ~Delaunay2D();

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
	void mesh(const TAlgoStrategy AStrategy = DELAUNAY_STRATEGY);

	/*------------------------------------------------------------------------*/
    /** \brief  Build a mesh on the points and segments stored in (*this)
	*/
	void meshNew(const TAlgoStrategy AStrategy = DELAUNAY_STRATEGY);

	/*------------------------------------------------------------------------*/
    /** \brief  return a reference on the resulting mesh
	*/
	gmds::IGMesh& getMesh() {return mesh_;}

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
	bool swap(gmds::Face& AF1, gmds::Face& AF2, const bool& AB = false);

	/*------------------------------------------------------------------------*/
    /** \brief  provide the strip of faces intersected by the segment [AN1,AN2]
     *
     *
	 */
	void getTriangleStrip(	gmds::Node& AN1, gmds::Node& AN2,
							std::vector<gmds::Face>& AStrip);

	/*------------------------------------------------------------------------*/
    /** \brief  constraint the triangular mesh to own the edge forming by nodes
     * 			AN1 and AN2
     *
     *
	 */
	void constraintToEdge(gmds::Node& AN1, gmds::Node& AN2);

	int getBoundaryMark() const {return init_mark_;}
protected:

	void middleEdgeInsertion();
	void middleEdgeInsertionOld();
	void edgeLengthBasedInsertion();
	/*------------------------------------------------------------------------*/
    /** \brief  class used to order 2D points in a std::set
	*/
	class PntComp{

	public:

		int operator()(const gmds::math::Point& p1,
		const gmds::math::Point& p2) const {
			return ((p1.X()<p2.X()) || ((p1.X()==p2.X()) && (p1.Y()<p2.Y())) );}
	};

	void addBoundary(const OrientedBoundary& ABoundary);
	/*------------------------------------------------------------------------*/
    /** \brief  Build a triangulation on the points stored in (*this)
	*/
	void initTriangulation(gmds::math::Point& AP_xmin_ymin,
			gmds::math::Point& AP_xmax_ymin,
			gmds::math::Point& AP_xmax_ymax,
			gmds::math::Point& AP_xmin_ymax);

	bool isInCircumCircle(gmds::Face& AFace,
						  gmds::math::Point& p);
	bool isInRelaxedCircumCircle(gmds::Face& AFace,
						  gmds::math::Point& p);

	void buildCavity(gmds::Face& AFace, gmds::math::Point& p,
					 std::vector<gmds::Face>& cavity,
					 const bool& relax=false);

	void splitCavity(gmds::Node& ANodeToInsert,
					 std::vector<gmds::Face>& cavity);
	void splitCavityOld(gmds::Node& ANodeToInsert,
					 std::vector<gmds::Face>& cavity);

	gmds::Face findFirstTriangle(gmds::math::Point& p);
	gmds::Face findRelaxedFirstTriangle(gmds::math::Point& p);
	void insertNode(gmds::Node& ANode);

	void removeOuterTrianglesNew();
	void removeOuterTriangles();
	void removeOuterTrianglesWithRays();
	/*------------------------------------------------------------------------*/
    /** \brief  Starting from an existing triangulation, this operation
     * 			constraint the triangulation to capture the edges of
     * 			init_segments_.
	*/
	void constraintBoundary();

	/*------------------------------------------------------------------------*/
    /** \brief  Check if the edge [AN1,AN2] is available in the triangulation
	*/
	bool checkEdgeExistence(gmds::Node& AN1, gmds::Node& AN2);

	/*------------------------------------------------------------------------*/
    /** \brief  Check if the point (x,y) is out of the domain
	*/
	bool isOutOfTheDomain(const TNum& AX, const TNum& AY);

	void displayMeshInfo();
	void initVariables();
	void initADimensionalSize();
	void setADimensionalSize(gmds::Face*);

	bool isBoundaryEdge(gmds::Node& AN1, gmds::Node& AN2);

	class compareFaceByRadius
	{
	 public:
	  inline bool operator () (const gmds::Face& f1, const gmds::Face& f2)  const
	  {
//		  if(adimensional_size_field_[f1->getID()]>adimensional_size_field_[f2->getID()])
//			  return true;
//		  else
			  return false;

	  }
	};

	void computeCircumCenter(gmds::Face& f, gmds::math::Point& center);

private:

	int init_mark_;

	std::vector<std::pair<gmds::Node, gmds::Node > > init_segments_;
	std::map<gmds::TCellID,std::vector<gmds::TCellID> > init_boundary_association_;
	std::vector<gmds::Node> bounding_box_nodes_;
	std::vector<gmds::Node> init_nodes_;
	std::vector<gmds::math::Point > init_points_;
	gmds::Variable<std::list<gmds::Face> >* init_N2F_node_variable;
	gmds::Variable<TNum >* size_field;

	/* variable associated to faces to store their adimensional size*/
	gmds::Variable<TNum >* adimensional_size_field_;
	gmds::Variable<gmds::math::Point >* circumCenter_field_;

	/* Boolean mark used to indicate if a face has been deleted or not from
	 * the mesh. */
	int delete_mark_;

	gmds::IGMesh mesh_;
};
/*----------------------------------------------------------------------------*/
#include "Delaunay2D.t.h"
/*----------------------------------------------------------------------------*/
}
/*----------------------------------------------------------------------------*/
#endif /* DELAUNAY2D_H_ */
/*----------------------------------------------------------------------------*/

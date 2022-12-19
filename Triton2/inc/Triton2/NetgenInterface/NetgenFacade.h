/*----------------------------------------------------------------------------*/
/*
 * NetgenFacade.h
 *
 *  Created on: 16 mai 2011
 *      Author: ledouxf
 */
/*----------------------------------------------------------------------------*/
#ifndef NetgenFACADE_H_
#define NetgenFACADE_H_
/*----------------------------------------------------------------------------*/
#include <GMDS/Utils/CommonTypes.h>
#include <GMDS/IG/IGMesh.h>
#include <GMDS/CAD/GeomVolume.h>
//#include <GMDSGeom/GeomManager.h>
#include <GMDS/CAD/GeomServiceAbstractFactory.h>
#include <GMDS/CAD/GeomTriangulationService.h>

/*----------------------------------------------------------------------------*/
#include <GMDS/Math/Point.h>
/*----------------------------------------------------------------------------*/
#include <sstream>
#include <map>
#include <set>
/*----------------------------------------------------------------------------*/
//inclusion de fichiers en-tête d'Open Cascade
#include <TopoDS_Edge.hxx>
#include <TDF_Label.hxx>
#include <gp_Pnt.hxx>
#include <GC_MakeSegment.hxx>
#include <TColgp_Array1OfPnt.hxx>
#include <BRepBuilderAPI_MakeEdge.hxx>
#include <ShapeConstruct.hxx>
#include <Geom_Curve.hxx>
#include <Geom_TrimmedCurve.hxx>
#include <TopoDS_Shape.hxx>
#include <TopoDS.hxx>
#include <ShapeExtend_WireData.hxx>
#include <ShapeFix_Wire.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <TopoDS_Face.hxx>
/*----------------------------------------------------------------------------*/
//inclusion de fichiers en-tête de Netgen
#include "NetgenLib/nglib/nglib.h"
#include "NetgenLib/nglib/nglib_addon.h"
//#include "nglib.h"
//#include "occgeom.hpp"
//#include "myadt.hpp"
using namespace nglib;
/*----------------------------------------------------------------------------*/
namespace triton{
/*----------------------------------------------------------------------------*/
/** \class NetgenParams
 *  \brief This class encapsulates the Netgen parameters to mesh a surface
 */
class NetgenParams{

public:

	bool uselocalh;              //!< Switch to enable / disable usage of local mesh size modifiers

	double maxh;                //!< Maximum global mesh size allowed
	double minh;                //!< Minimum global mesh size allowed

	double fineness;            //!< Mesh density: 0...1 (0 => coarse; 1 => fine)
	double grading;             //!< Mesh grading: 0...1 (0 => uniform mesh; 1 => aggressive local grading)

	double elementsperedge;     //!< Number of elements to generate per edge of the geometry
	double elementspercurve;    //!< Elements to generate per curvature radius

	bool closeedgeenable;        //!< Enable / Disable mesh refinement at close edges
	double closeedgefact;       //!< Factor to use for refinement at close edges
								//!<(STL: larger => finer ; OCC: larger => coarser)

	bool secondorder;            //!< Generate second-order surface and volume elements
	bool quad_dominated;         //!< Creates a Quad-dominated mesh

	bool optsurfmeshenable;      //!< Enable / Disable automatic surface mesh optimization
	bool optvolmeshenable;       //!< Enable / Disable automatic volume mesh optimization

	int optsteps_3d;            //!< Number of optimize steps to use for 3-D mesh optimization
	int optsteps_2d;            //!< Number of optimize steps to use for 2-D mesh optimization


	/*!
	      Default constructor for the Mesh Parameters class
	 */
	NetgenParams(){
		uselocalh=true;
		maxh= 1000.0;
		minh = 0.01;
		fineness=0.5;
		grading = 0.3;
		elementsperedge= 10.0;
	    elementspercurve=10.0;
	    closeedgeenable= false;
	    closeedgefact= 2.0;
	    secondorder= false;
	    quad_dominated= false;
	    optsurfmeshenable=true;
	    optvolmeshenable=true;
	    optsteps_2d= 3;
	    optsteps_3d= 3;
	}
};
/*----------------------------------------------------------------------------*/
/** \class NetgenFacade
 *  \brief This class encapsulates the different calls to Netgen.
 */
class NetgenFacade{

public:
	/*------------------------------------------------------------------------*/
    /** \brief  Constructor.
     *
     *	\param ABoundary the outter boundary of the domain we want to mesh.
     */
	NetgenFacade();


	/*------------------------------------------------------------------------*/
    /** \brief  This method allows to tetrahedralize a volume defined by the
     * 			boundary nodes and faces provided as arguments. These nodes and
     * 			faces have to belong to AVolMesh which is the mesh in which
     * 			the resulting tetrahedral elements will be added.
     *
     *
     *	\param ANodes		boundary nodes
     *	\param AFaces		boundary faces
     *	\param AVolMesh		the resulting volume mesh
     */
	void generateTetMesh(std::vector<gmds::Node>& ANodes,
			std::vector<gmds::Face>& AFaces,
			gmds::IGMesh& AVolMesh);

	void generateTetMesh(std::vector<gmds::Node>& ANodes,
			std::vector<gmds::Face>& AFaces,
			std::vector<gmds::Node>& ACreatedNodes,
			std::vector<gmds::Region>& ACreatedRegions,
			gmds::IGMesh& AVolMesh,
			NetgenParams& AParams);

	void generateTetMesh(gmds::IGMesh& AVolMesh, const double ASize);

	/*------------------------------------------------------------------------*/
    /** \brief  This method allows to tetrahedralize a volume defined by the
     * 			boundary mesh ABoundaryMesh. The result will be added in
     * 			AVolMesh. In order to get the expected AQuality, boundary nodes
     * 		    can be added.The radius-edge ratio must be included in [1.0,
     * 			2.0].
     *
     *
     *	\param ABoundaryMesh the outter boundary of the domain we want to mesh.
     *	\param AVolMesh		 the resulting volume mesh
     *	\param AQuality 	 defines a quality target size for each tetrahedron
     */
	void generateTetMesh(gmds::IGMesh& ABoundaryMesh,
			gmds::IGMesh& AVolMesh,
						const double AQuality = 2.0);

	/*------------------------------------------------------------------------*/
    /** \brief  This method allows us to triangulate a volume defined by an
     * 			object implementing the gmds::geom::GeomVolume interface.
     * 			The result will be added in ASurfMesh. In order to get the
     * 			expected AQuality, boundary nodes can be added.The radius-edge
     * 			ratio must be included in [1.0,	2.0].
     *
     *
     *	\param AVol			the geometric volume we want to mesh
     *	\param ASurfMesh		the resulting surface mesh
     *	\param AQuality 	defines a quality target size for each tetrahedron
     */
	void generateTriMesh(gmds::geom::GeomVolume& AVol,
						 gmds::geom::GeomServiceAbstractFactory& AMan,
						 gmds::IGMesh& ASurMesh,
						 const double AQuality = 2.0);

	void generateTriMesh(gmds::geom::GeomVolume& AVol,
						 gmds::geom::GeomServiceAbstractFactory& AMan,
						 gmds::IGMesh& ASurMesh,
						 const double AQuality,
						 const double ASize);

	/*------------------------------------------------------------------------*/
    /** \brief  This method allows us to triangulate a volume defined by an
     * 			object implementing the gmds::geom::GeomVolume interface.
     * 			The result will be added in ASurfMesh. In order to get the
     * 			expected AQuality, boundary nodes can be added.The radius-edge
     * 			ratio must be included in [1.0,	2.0].
     *
     *
     *	\param ASurf		the geometric surface we want to mesh
     *	\param ASurfMesh	the resulting surface mesh
     *	\param AQuality 	defines a quality target size for each tetrahedron
     */
	void generateTriMesh(gmds::geom::GeomSurface& ASurf,
						 gmds::geom::GeomServiceAbstractFactory& AMan,
						 gmds::IGMesh& ASurMesh,
						 const double AQuality = 2.0);

	/*------------------------------------------------------------------------*/
    /** \brief  This method allows us to tetrahedralize a volume defined by an
     * 			object implementing the gmds::geom::GeomVolume interface.
     * 			The result will be added in AVolMesh. In order to get the
     * 			expected AQuality, boundary nodes can be added.The radius-edge
     * 			ratio must be included in [1.0,	2.0].
     *
     *
     *	\param AVol			the geometric volume we want to mesh
     *	\param AVolMesh		the resulting volume mesh
     *	\param AQuality 	defines a quality target size for each tetrahedron
     */
	void generateTetMesh(gmds::geom::GeomVolume& AVol,
						 gmds::geom::GeomServiceAbstractFactory& AMan,
						 gmds::IGMesh& AVolMesh,
						 const double AQuality = 2.0);
	void generateTetMesh(gmds::geom::GeomVolume& AVol,
						 gmds::geom::GeomServiceAbstractFactory& AMan,
						 gmds::IGMesh& AVolMesh,
						 const double AQuality,
						 gmds::math::Point& APMin,
						 gmds::math::Point& APMax,
						 const double ASize);
	/*------------------------------------------------------------------------*/
    /** \brief  Initializes the Netgen internal data structure from AMesh, a
     * 			GMDS mesh given in parameter
     *
     *	\param AMesh the input mesh
     *	\param ASTLGeom the STL Geometry in the Netgen API
     */
	void generateTriMeshFromSTEP(const std::string& AName,
			 	 	 	 	 	 gmds::IGMesh& ASurfMesh,
								 const double AQuality = 2.0);

	/*------------------------------------------------------------------------*/
    /** \brief  Initializes the Netgen internal data structure from AMesh, a
     * 			GMDS mesh given in parameter
     *
     *	\param AMesh the input mesh
     *	\param ASTLGeom the STL Geometry in the Netgen API
     */
	void generateTriMesh(const std::string& AName,
						 gmds::IGMesh& ASurfMesh,
						 std::vector<gmds::Node>& ANodes,
						 std::vector<std::vector<gmds::Node> >& AEdges,
						 const double AQuality = 2.0);

	void generateTriMesh(TopoDS_Face& AFace,
						 gmds::IGMesh& ASurfMesh,
						 std::vector<gmds::Node>& ANodes,
						 std::vector<std::vector<gmds::Node> >& AEdges,
						 std::map<gmds::TCellID,gmds::TCellID>& ANodeMap,
						  NetgenParams& AParams);
	/*------------------------------------------------------------------------*/
    void generateTriMeshWithOCC( std::vector<gmds::math::Point >& APnts,
				 	 	 	 	 gmds::IGMesh& ASurfMesh,
								 const double AQuality = 2.0);

	TopoDS_Face buildOCCSurfaceFromPoints(std::vector<gmds::math::Point >& APnts);
	void meshing2D( TopoDS_Face & AFace);
private:

	/*------------------------------------------------------------------------*/
    /** \brief  Initializes the Netgen parameters and create Netgen objects
     */
	void initializeNetgen(NetgenParams& AParams);
	void initializeNetgen(const double AQuality,
			gmds::math::Point& APMin,
			gmds::math::Point& APMax,
			 const double ASize);
	void initializeNetgen(const double AQuality,
			 const double ASize=1);
	/*------------------------------------------------------------------------*/
    /** \brief  Clean the Netgen objects
     */
	void finalizeNetgen();

	/*------------------------------------------------------------------------*/
    /** \brief  Initializes the Netgen internal data structure from AMesh, a
     * 			GMDS mesh given in parameter
     *
     *	\param AMesh the input mesh
     *	\param ASTLGeom the STL Geometry in the Netgen API
     */
	void buildSTLSurfaceFrom(gmds::IGMesh& AMesh);

	/*------------------------------------------------------------------------*/
    /** \brief  Initializes the Netgen internal data structure from AVol, a
     * 			geometric volume we want to mesh
     *
     *	\param AVol	 	a geometric volume
     *	\param ASTLGeom the STL Geometry in the Netgen API
     */
	void buildSTLSurfaceFrom(gmds::geom::GeomVolume& AVol,
			 gmds::geom::GeomServiceAbstractFactory& AMan);
	/*------------------------------------------------------------------------*/
    /** \brief  Initializes the Netgen internal data structure from AVol, a
     * 			geometric volume we want to mesh
     *
     *	\param ASurf	 	a geometric surface
     *	\param ASTLGeom the STL Geometry in the Netgen API
     */
	void buildSTLSurfaceFrom(gmds::geom::GeomSurface& ASurf,
			 gmds::geom::GeomServiceAbstractFactory& AMan);

	/*------------------------------------------------------------------------*/
    /** \brief  Atomic operation to create the expected volume mesh
     *
     *	\param AMesh the output mesh
     */
	void createGMDSMesh(gmds::IGMesh& AMesh);
	void createGMDSMesh(gmds::IGMesh& AMesh,
			gmds::math::Point& APMin,
			gmds::math::Point& APMax,
			 const double ASize);

	/*------------------------------------------------------------------------*/
    /** \brief  Atomic operation to create the expected volume mesh
     *
     *	\param AMesh the output mesh
     */
	void createGMDSTriMesh(gmds::IGMesh& AMesh);

	/*------------------------------------------------------------------------*/
    /** \brief  Build AMesh from the NetgenOutput_ data structure
     *
     *	\param AMesh the input mesh
     */
	void buildGMDSOutput(gmds::IGMesh& AMesh);
	void buildGMDSOutput(gmds::IGMesh& AMesh,std::set<gmds::TCellID>& initNodes,
						std::map<gmds::TCellID,gmds::TCellID>& ANodeMap);

	Ng_Mesh* ng_mesh_;
	Ng_Meshing_Parameters ng_mesh_param_surf_;
	Ng_Meshing_Parameters ng_mesh_param_vol_;
	Ng_STL_Geometry* stl_geometry_;
};

/*----------------------------------------------------------------------------*/
}
/*----------------------------------------------------------------------------*/
#endif /* NETGENFACADE_H_ */
/*----------------------------------------------------------------------------*/

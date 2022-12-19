/*----------------------------------------------------------------------------*/
/*
 * GeomToolKit.h
 *
 *  Created on: 19 janv. 2011
 *      Author: ledouxf
 */
/*----------------------------------------------------------------------------*/
#ifndef GEOMTOOLKIT_H_
#define GEOMTOOLKIT_H_
/*----------------------------------------------------------------------------*/
#include "GMDS/Math/Point.h"
#include "GMDS/Math/Matrix.h"
/*--------------------------------------------------------------------------*/
#include "GMDS/Utils/CommonTypes.h"
/*--------------------------------------------------------------------------*/
#include <vector>
/*----------------------------------------------------------------------------*/
namespace triton{
/*----------------------------------------------------------------------------*/
/** \class GeomToolKit
 *  \brief Gathers geometric algorithms necessary to perform Delaunay
 *  	   triangulation. In this context, mesh notions are unknown.
 */
template <typename T>
class GeomToolKit{

public:

	/*------------------------------------------------------------------------*/
    /** \brief  Predicate if a and b are equal with respect to tolerance
	*/
	static bool compare(const T& a, const T& b,
					    const T& tolerance = 1e-12);
	/*------------------------------------------------------------------------*/
    /** \brief  Compute the bounding box of APnts with tolerance ATol
	*/
	static void getBoundingBox(
			const std::vector<gmds::math::Point >& APnts,
			gmds::math::Point& APmin,
			gmds::math::Point& APmax,
			const T& ATol = 0.1);

	/*------------------------------------------------------------------------*/
    /** \brief  Compute the circumcenter of 3 2D points
	*/
	static void computeCircumCenter(
			const gmds::math::Point& p1, const gmds::math::Point& p2,
			const gmds::math::Point& p3, gmds::math::Point& center);

	/*------------------------------------------------------------------------*/
    /** \brief  Predicate indication if point n is in the circumcircle of the
     * 			triangle defining by p1, p2, p3. Points p1, p2 and p3 have not
     * 			to be counterclockwise ordered.
     *
     * 			It returns a positive value if n lies inside the circle passing
     * 			through p1, p2, and p3; a negative value if it lies outside;
     * 			and zero if the five points are cocircular.
	*/
	static T isInCircumCircle(
			gmds::math::Point& p1, gmds::math::Point& p2,
			gmds::math::Point& p3, gmds::math::Point& n);

	/*------------------------------------------------------------------------*/
    /** \brief  Computes the shortest distance from a point AP to a line
     * 			defined by points AP1 and AP2
     	*/
	static T shortestDistance(gmds::math::Point& AP,
			gmds::math::Point& AP1,
			gmds::math::Point& AP2);
};
/*----------------------------------------------------------------------------*/
#include "GeomToolKit.t.h"
/*----------------------------------------------------------------------------*/
}
/*----------------------------------------------------------------------------*/
#endif /* GEOMTOOLKIT_H_ */
/*----------------------------------------------------------------------------*/


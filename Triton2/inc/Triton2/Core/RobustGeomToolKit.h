/*----------------------------------------------------------------------------*/
/*
 * RobustGeomToolKit.h
 *
 *  Created on: 21 janv. 2011
 *      Author: ledouxf
 */
/*----------------------------------------------------------------------------*/
#ifndef ROBUSTGEOMTOOLKIT_H_
#define ROBUSTGEOMTOOLKIT_H_
/*----------------------------------------------------------------------------*/
#include <GMDS/Utils/CommonTypes.h>
#include <GMDS/Math/Point.h>
#include <vector>
/*----------------------------------------------------------------------------*/
namespace triton{
/*----------------------------------------------------------------------------*/
/** \class RobustGeomToolKit
 *  \brief Gathers geometric algorithms necessary to perform Delaunay
 *  	   triangulation. In this context, mesh notions are unknown.
 *
 *  	   	The implementation of the geometric predicates are those of J.
 *  		Shewchuk
 */
class RobustGeomToolKit{

public:

	/*------------------------------------------------------------------------*/
    /** \brief  Predicate if a and b are equal with respect to tolerance
	*/
	static bool compare(const gmds::TCoord& a, const gmds::TCoord& b,
					    const gmds::TCoord& tolerance = 1e-12);
	/*------------------------------------------------------------------------*/
    /** \brief  Compute the bounding box of APnts with tolerance ATol
	*/
	static void getBoundingBox(
			const std::vector<gmds::math::Point >& APnts,
			gmds::math::Point& APmin,
			gmds::math::Point& APmax,
			const gmds::TCoord& ATol = 0.1);
	/*------------------------------------------------------------------------*/
    /** \brief  Predicate indication if point n is in the circumcircle of the
     * 			triangle defining by p1, p2, p3. Points p1, p2 and p3 have not
     * 			to be counterclockwise ordered.
     *
     * 			It returns a positive value if n lies inside the circle passing
     * 			through p1, p2, and p3; a negative value if it lies outside;
     * 			and zero if the five points are cocircular.
	*/
	static gmds::TCoord isInCircumCircle(
			gmds::math::Point& p1, gmds::math::Point& p2,
			gmds::math::Point& p3, gmds::math::Point& n);

};
/*----------------------------------------------------------------------------*/
}
/*----------------------------------------------------------------------------*/
#endif /* ROBUSTGEOMTOOLKIT_H_ */
/*----------------------------------------------------------------------------*/



/*----------------------------------------------------------------------------*/
/*
 * robustPredicates.h
 *
 *  Created on: 21 janv. 2011
 *      Author: ledouxf
 */
/*----------------------------------------------------------------------------*/
#ifndef ROBUSTPREDICATES_H_
#define ROBUSTPREDICATES_H_
/*----------------------------------------------------------------------------*/
namespace triton{
/*----------------------------------------------------------------------------*/
/** \class RobustPredicates
 *  \brief This class is a wrapper to the robust predicates defined and
 *  	   implemented by J. Shewchuk for Triangle.
 */
class RobustPredicates{

public:

	/*------------------------------------------------------------------------*/
    /** \brief  Initialize the variables used for exact arithmetic.
	*/
	static double exactinit();
	/*------------------------------------------------------------------------*/
    /** \brief  Adaptive exact 2D incircle test.
	/*
	/*          Returns a positive value if the point pd lies inside the
	/*          circle passing through pa, pb, and pc; a negative value if
	/*          it lies outside; and zero if the four points are cocircular.
	/*          The points pa, pb, and pc must be in counterclockwise
	/*          order, or the sign of the result will be reversed.
	*/
	static double incircle(double *pa, double *pb, double *pc, double *pd);
	/*------------------------------------------------------------------------*/
    /** \brief  Adaptive exact 3D insphere test.
	/*
	/*          Returns a positive value if the point pe lies inside the
	/*          sphere passing through pa, pb, pc, and pd; a negative value
	/*          if it lies outside; and zero if the five points are
	/*          cospherical.  The points pa, pb, pc, and pd must be ordered
	/*          so that they have a positive orientation (as defined by
	/*          orient3d()), or the sign of the result will be reversed.
	*/
	static double insphere(double *pa, double *pb, double *pc, double *pd,
						   double *pe);
	/*------------------------------------------------------------------------*/
    /** \brief  Adaptive exact 2D orientation test.
	/*
	/*          Returns a positive value if the points pa, pb, and pc occur
	/*          in counterclockwise order; a negative value if they occur
	/*          in clockwise order; and zero if they are collinear.  The
	/*          result is also a rough approximation of twice the signed
	/*          area of the triangle defined by the three points.
	*/
	static double orient2d(double *pa, double *pb, double *pc);
	/*------------------------------------------------------------------------*/
    /** \brief  Adaptive exact 3D orientation test.
	/*
	/*          Returns a positive value if the point pd lies below the
	/*          plane passing through pa, pb, and pc; "below" is defined so
	/*          that pa, pb, and pc appear in counterclockwise order when
	/*          viewed from above the plane.  Returns a negative value if
	/*          pd lies above the plane.  Returns zero if the points are
	/*          coplanar.  The result is also a rough approximation of six
	/*          times the signed volume of the tetrahedron defined by the
	/*          four points.
	*/
	static double orient3d(double *pa, double *pb, double *pc, double *pd);
};
/*----------------------------------------------------------------------------*/
}
/*----------------------------------------------------------------------------*/
#endif /* ROBUSTPREDICATES_H_ */
/*----------------------------------------------------------------------------*/



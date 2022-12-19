/*----------------------------------------------------------------------------*/
/*
 * GeomToolKitTestSuite.cpp
 *
 *  Created on: 19 janv. 2011
 *      Author: ledouxf
 */
/*----------------------------------------------------------------------------*/
#include "GeomToolKitTestSuite.h"
#include "Triton2/Core/GeomToolKit.h"
#include <math.h>
/*----------------------------------------------------------------------------*/
using namespace gmds;
using namespace triton;
/*----------------------------------------------------------------------------*/
void GeomToolKitTestSuite::testInCircumCircle1(){
	math::Point p1(0.0,0.0), p2(1.0,0.0), p3(0.0,1.0), n(1.0,1.0);

	gmds::TCoord r =GeomToolKit<gmds::TCoord>::isInCircumCircle(p1,p2,p3,n);
	CPPUNIT_ASSERT(r==0.0);
}
/*----------------------------------------------------------------------------*/
void GeomToolKitTestSuite::testInCircumCircle2(){
	math::Point p1(0.0,0.0), p2(1.0,0.0), p3(0.0,1.0), n(0.25,0.25);
	gmds::TCoord r =GeomToolKit<gmds::TCoord>::isInCircumCircle(p1,p2,p3,n);
	CPPUNIT_ASSERT(r>0.0);
}
/*----------------------------------------------------------------------------*/
void GeomToolKitTestSuite::testInCircumCircle3(){
	math::Point p1(0.0,0.0), p2(1.0,0.0), p3(0.0,1.0), n(0.5,0.5);
	gmds::TCoord r =GeomToolKit<gmds::TCoord>::isInCircumCircle(p1,p2,p3,n);
	CPPUNIT_ASSERT(r>0.0);
}
/*----------------------------------------------------------------------------*/
void GeomToolKitTestSuite::testInCircumCircle4(){
	math::Point p1(0.0,0.0), p2(1.0,0.0), p3(0.0,1.0), n(0.5,0.6);
	gmds::TCoord r =GeomToolKit<gmds::TCoord>::isInCircumCircle(p1,p2,p3,n);
    CPPUNIT_ASSERT(r>0.0);
}
/*----------------------------------------------------------------------------*/
void GeomToolKitTestSuite::testInCircumCircle5(){
	math::Point p1(0.0,0.0), p2(1.0,0.0), p3(0.0,1.0), n(0.0,-0.01);
	gmds::TCoord r =GeomToolKit<gmds::TCoord>::isInCircumCircle(p1,p2,p3,n);
	CPPUNIT_ASSERT(r<0.0);
}
/*----------------------------------------------------------------------------*/
void GeomToolKitTestSuite::testInCircumCircle6(){
	math::Point p1(0.0,0.0), p2(1.0,0.0), p3(0.0,1.0);
	math::Point c(p1.X()+p2.X()+p3.X()/3, p1.Y()+p2.Y()+p3.Y()/3);
	gmds::TCoord radius = sqrt((p1.X()-c.X())*(p1.X()-c.X())+ (p1.Y()-c.Y())*(p1.Y()-c.Y()));

	math::Point n(c.X()+sqrt(radius),c.Y()+sqrt(radius));
	gmds::TCoord r =GeomToolKit<gmds::TCoord>::isInCircumCircle(p1,p2,p3,n);
	CPPUNIT_ASSERT(r<0.0);
}
/*----------------------------------------------------------------------------*/
void GeomToolKitTestSuite::testInCircumCircle7(){
	math::Point p1(0.0,0.0), p2(1.0,0.0), p3(0.0,1.0);
	gmds::TCoord r =GeomToolKit<gmds::TCoord>::isInCircumCircle(p1,p2,p3,p1);
	CPPUNIT_ASSERT(r==0.0);
}
/*----------------------------------------------------------------------------*/
void GeomToolKitTestSuite::testInCircumCircle8(){
	math::Point p1(0.0,0.0), p2(1.0,0.0), p3(0.0,1.0);
	gmds::TCoord r =GeomToolKit<gmds::TCoord>::isInCircumCircle(p1,p2,p3,p2);
	CPPUNIT_ASSERT(r==0.0);
}
/*----------------------------------------------------------------------------*/
void GeomToolKitTestSuite::testInCircumCircle9(){
	math::Point p1(0.0,0.0), p2(1.0,0.0), p3(0.0,1.0);
	gmds::TCoord r =GeomToolKit<gmds::TCoord>::isInCircumCircle(p1,p2,p3,p3);
	CPPUNIT_ASSERT(r==0.0);
}
/*----------------------------------------------------------------------------*/
void GeomToolKitTestSuite::testInCircumCircle10(){
	math::Point p1(0.456,2.546), p2(-3.41567,0.78945612),
						  p3(4.014567,-0.7896247);
	gmds::TCoord r =GeomToolKit<gmds::TCoord>::isInCircumCircle(p1,p2,p3,p1);
	CPPUNIT_ASSERT(r==0.0);
}
/*----------------------------------------------------------------------------*/
void GeomToolKitTestSuite::testInCircumCircle11(){
	math::Point p1(0.456,2.546), p2(-3.41567,0.78945612),
						  p3(4.014567,-0.7896247);
	gmds::TCoord r =GeomToolKit<gmds::TCoord>::isInCircumCircle(p1,p2,p3,p2);
	CPPUNIT_ASSERT(GeomToolKit<gmds::TCoord>::compare(r,0.0));
}
/*----------------------------------------------------------------------------*/
void GeomToolKitTestSuite::testInCircumCircle12(){
	math::Point p1(0.456,2.546), p2(-3.41567,0.78945612),
						  p3(4.014567,-0.7896247);
	gmds::TCoord r =GeomToolKit<gmds::TCoord>::isInCircumCircle(p1,p2,p3,p3);
	CPPUNIT_ASSERT(r==0.0);
}
/*----------------------------------------------------------------------------*/

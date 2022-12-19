/*----------------------------------------------------------------------------*/
/*
 * RobustRobustGeomToolKitTestSuite.cpp
 *
 *  Created on: 19 janv. 2011
 *      Author: ledouxfgmds::TCoord
 */
/*----------------------------------------------------------------------------*/
#include "RobustGeomToolKitTestSuite.h"
#include "Triton2/Core/RobustGeomToolKit.h"
#include <math.h>
/*----------------------------------------------------------------------------*/
using namespace gmds;
using namespace triton;
/*----------------------------------------------------------------------------*/
void RobustGeomToolKitTestSuite::testInCircumCircle1(){
//	GEPETO::Point<2,gmds::TCoord> p1(0.0,0.0), p2(1.0,0.0), p3(1.0,1.0), n(1,1);
//	gmds::TCoord r = RobustGeomToolKit::isInCircumCircle(p1,p2,p3,n);
//	CPPUNIT_ASSERT(r==0.0);
}
/*----------------------------------------------------------------------------*/
void RobustGeomToolKitTestSuite::testInCircumCircle2(){
//	GEPETO::Point<2,gmds::TCoord> p1(0.0,0.0), p2(1.0,0.0), p3(1.0,1.0), n(0.25,0.25);
//	gmds::TCoord r =RobustGeomToolKit::isInCircumCircle(p1,p2,p3,n);
//	CPPUNIT_ASSERT(r>0.0);
}
/*----------------------------------------------------------------------------*/
void RobustGeomToolKitTestSuite::testInCircumCircle3(){
//	GEPETO::Point<2,gmds::TCoord> p1(0.0,0.0), p2(1.0,0.0), p3(1.0,1.0), n(0.5,0.5);
//	gmds::TCoord r =RobustGeomToolKit::isInCircumCircle(p1,p2,p3,n);
	//	CPPUNIT_ASSERT(r>0.0);
}
/*----------------------------------------------------------------------------*/
void RobustGeomToolKitTestSuite::testInCircumCircle4(){
//	GEPETO::Point<2,gmds::TCoord> p1(0.0,0.0), p2(1.0,0.0), p3(1.0,1.0), n(0.5,0.6);
//	gmds::TCoord r =RobustGeomToolKit::isInCircumCircle(p1,p2,p3,n);
	//CPPUNIT_ASSERT(r>0.0);
}
/*----------------------------------------------------------------------------*/
void RobustGeomToolKitTestSuite::testInCircumCircle5(){
//	GEPETO::Point<2,gmds::TCoord> p1(0.0,0.0), p2(1.0,0.0), p3(0.1,1.0), n(0,-0.01);
//	gmds::TCoord r =RobustGeomToolKit::isInCircumCircle(p1,p2,p3,n);
	//CPPUNIT_ASSERT(r<0.0);
}
/*----------------------------------------------------------------------------*/
void RobustGeomToolKitTestSuite::testInCircumCircle6(){
//	GEPETO::Point<2,gmds::TCoord> p1(0.0,0.0), p2(1.0,0.0), p3(1.0,1.0);
//	GEPETO::Point<2,gmds::TCoord> c(p1.getX()+p2.getX()+p3.getX()/3,
//							p1.getY()+p2.getY()+p3.getY()/3);
//	gmds::TCoord radius = gmds::TCoord::sqrt((p1.getX()-c.getX())*(p1.getX()-c.getX())+
//						 (p1.getY()-c.getY())*(p1.getY()-c.getY()));
//
//	GEPETO::Point<2,gmds::TCoord> n(c.getX()+gmds::TCoord::sqrt(radius),c.getY()+gmds::TCoord::sqrt(radius));
//	gmds::TCoord r =RobustGeomToolKit::isInCircumCircle(p1,p2,p3,n);
	//	CPPUNIT_ASSERT(r<0.0);
}
/*----------------------------------------------------------------------------*/
void RobustGeomToolKitTestSuite::testInCircumCircle7(){
//	GEPETO::Point<2,gmds::TCoord> p1(0.0,0.0), p2(1.0,0.0), p3(1.0,1.0);
//	gmds::TCoord r =RobustGeomToolKit::isInCircumCircle(p1,p2,p3,p1);
	//	CPPUNIT_ASSERT(r==0.0);
}
/*----------------------------------------------------------------------------*/
void RobustGeomToolKitTestSuite::testInCircumCircle8(){
//	GEPETO::Point<2,gmds::TCoord> p1(0.0,0.0), p2(1.0,0.0), p3(1.0,1.0);
//	gmds::TCoord r =RobustGeomToolKit::isInCircumCircle(p1,p2,p3,p2);
	//	CPPUNIT_ASSERT(r==0.0);
}
/*----------------------------------------------------------------------------*/
void RobustGeomToolKitTestSuite::testInCircumCircle9(){
//	GEPETO::Point<2,gmds::TCoord> p1(0.0,0.0), p2(1.0,0.0), p3(1.0,1.0);
//	gmds::TCoord r =RobustGeomToolKit::isInCircumCircle(p1,p2,p3,p3);
	//	CPPUNIT_ASSERT(r==0.0);
}
/*----------------------------------------------------------------------------*/
void RobustGeomToolKitTestSuite::testInCircumCircle10(){
//	GEPETO::Point<2,gmds::TCoord> p1(0.456,2.546), p2(-3.41567,0.78945612),
//						  p3(4.014567,-0.7896247);
//	gmds::TCoord r =RobustGeomToolKit::isInCircumCircle(p1,p2,p3,p1);
	//	CPPUNIT_ASSERT(r==0.0);
}
/*----------------------------------------------------------------------------*/
void RobustGeomToolKitTestSuite::testInCircumCircle11(){
//	GEPETO::Point<2,gmds::TCoord> p1(0.456,2.546), p2(-3.41567,0.78945612),
//						  p3(4.014567,-0.7896247);
//	gmds::TCoord r =RobustGeomToolKit::isInCircumCircle(p1,p2,p3,p2);
	//	CPPUNIT_ASSERT(RobustGeomToolKit::compare(r,0));
}
/*----------------------------------------------------------------------------*/
void RobustGeomToolKitTestSuite::testInCircumCircle12(){
//	GEPETO::Point<2,gmds::TCoord> p1(0.456,2.546), p2(-3.41567,0.78945612),
//						  p3(4.014567,-0.7896247);
//	gmds::TCoord r =RobustGeomToolKit::isInCircumCircle(p1,p2,p3,p3);
	//CPPUNIT_ASSERT(r==0.0);
}
/*----------------------------------------------------------------------------*/

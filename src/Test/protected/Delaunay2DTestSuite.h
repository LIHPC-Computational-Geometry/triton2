/*----------------------------------------------------------------------------*/
/*
 * Delaunay2DTestSuite.h
 *
 *  Created on: 19 janv. 2011
 *      Author: ledouxf
 */
/*----------------------------------------------------------------------------*/
#ifndef DELAUNAY2DTESTSUITE_H_
#define DELAUNAY2DTESTSUITE_H_
/*----------------------------------------------------------------------------*/
#include <string>
#include <cppunit/extensions/HelperMacros.h>
/*----------------------------------------------------------------------------*/
class Delaunay2DTestSuite: public CppUnit::TestFixture {

    CPPUNIT_TEST_SUITE(Delaunay2DTestSuite);
//    CPPUNIT_TEST(testInitDelaunay);
//    CPPUNIT_TEST(testDelaunaySquare);
//    CPPUNIT_TEST(testDelaunayUnitSquare);
//    CPPUNIT_TEST(testDelaunayEmptySquare);
   // CPPUNIT_TEST(testDelaunayCircle);
//    CPPUNIT_TEST(testNonConvex);
//    CPPUNIT_TEST(testNonConvex2);
//    CPPUNIT_TEST(testNonConvex3);
//    CPPUNIT_TEST(testNonConvex4);
    CPPUNIT_TEST(testDoubleCircles);
//    CPPUNIT_TEST(testSwapConvex);
//    CPPUNIT_TEST(testSwapNonConvex);
//    CPPUNIT_TEST(testGetTriangleStrip);
//    CPPUNIT_TEST(testBuildEdgeWith6ConvexTriangles);
//    CPPUNIT_TEST(testBuildEdgeWith3NonConvexTriangles);
    CPPUNIT_TEST_SUITE_END();

public:
    void setUp(){}

    void tearDown(){}

    void testInitDelaunay();

    void testDelaunaySquare();

    void testDelaunayUnitSquare();

    void testDelaunayEmptySquare();

    void testDelaunayCircle();

    void testNonConvex();
    void testNonConvex2();
    void testNonConvex3();
    void testNonConvex4();
    void testDoubleCircles();

    void testSwapConvex();
    void testSwapNonConvex();
    void testGetTriangleStrip();
    void testBuildEdgeWith6ConvexTriangles();
//    void testBuildEdgeWith3NonConvexTriangles();
};
/*----------------------------------------------------------------------------*/
#endif /* DELAUNAY2DTESTSUITE_H_ */
/*----------------------------------------------------------------------------*/



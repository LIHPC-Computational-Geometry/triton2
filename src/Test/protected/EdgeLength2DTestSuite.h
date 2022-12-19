/*----------------------------------------------------------------------------*/
/*
 * EdgeLength2DTestSuite.h
 *
 *  Created on: 29 oct. 2012
 *      Author: ledouxf
 */
/*----------------------------------------------------------------------------*/
#ifndef EDGELENGTH2DTESTSUITE_H_
#define EDGELENGTH2DTESTSUITE_H_
/*----------------------------------------------------------------------------*/
#include <cppunit/extensions/HelperMacros.h>
/*----------------------------------------------------------------------------*/
class EdgeLength2DTestSuite: public CppUnit::TestFixture {

    CPPUNIT_TEST_SUITE(EdgeLength2DTestSuite);
    CPPUNIT_TEST(testDoubleCircles);
    CPPUNIT_TEST_SUITE_END();

public:
    void setUp(){}
    void tearDown(){}
    void testDoubleCircles();
};
/*----------------------------------------------------------------------------*/
#endif /* EDGELENGTH2DTESTSUITE_H_ */
/*----------------------------------------------------------------------------*/





/*----------------------------------------------------------------------------*/
/*
 * NetgenFacadeTestSuite.h
 *
 *  Created on: 23 mai 2011
 *      Author: ledouxf
 */
/*----------------------------------------------------------------------------*/
#ifndef NETGENFACADETESTSUITE_H_
#define NETGENFACADETESTSUITE_H_
/*----------------------------------------------------------------------------*/
#include <string>
#include <cppunit/extensions/HelperMacros.h>
/*----------------------------------------------------------------------------*/
class NetgenFacadeTestSuite: public CppUnit::TestFixture {

    CPPUNIT_TEST_SUITE(NetgenFacadeTestSuite);
    CPPUNIT_TEST(testTet);
    CPPUNIT_TEST(testCube);
    CPPUNIT_TEST(testSurface);
    CPPUNIT_TEST(testSTEPSurface);
    CPPUNIT_TEST(testOCCFaceAndMeshNodes);
    CPPUNIT_TEST(testSTEPSolid);
    CPPUNIT_TEST(testOCCSurface);
    CPPUNIT_TEST(testM3DBlock);
    CPPUNIT_TEST(testFacetModel);
    CPPUNIT_TEST_SUITE_END();

public:

    void testTet();
    void testCube();
    void testSurface();
    void testSTEPSolid();
    void testSTEPSurface();
    void testM3DBlock();
    void testOCCSurface();
    void testFacetModel();
    void testOCCFaceAndMeshNodes();
};
/*----------------------------------------------------------------------------*/
#endif /* NETGENFACADETESTSUITE_H_ */
/*----------------------------------------------------------------------------*/


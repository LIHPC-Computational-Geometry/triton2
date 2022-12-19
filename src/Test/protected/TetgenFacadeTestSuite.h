/*----------------------------------------------------------------------------*/
/*
 * TetgenFacadeTestSuite.h
 *
 *  Created on: 17 mai 2011
 *      Author: ledouxf
 */
/*----------------------------------------------------------------------------*/
#ifndef TETGENFACADETESTSUITE_H_
#define TETGENFACADETESTSUITE_H_
/*----------------------------------------------------------------------------*/
#include <string>
#include <cppunit/extensions/HelperMacros.h>
/*----------------------------------------------------------------------------*/
class TetgenFacadeTestSuite: public CppUnit::TestFixture {

    CPPUNIT_TEST_SUITE(TetgenFacadeTestSuite);
    CPPUNIT_TEST(testM3DBlock);
//    CPPUNIT_TEST(testPLCCube);
//    CPPUNIT_TEST(testPLCCube2);
    CPPUNIT_TEST_SUITE_END();

public:
    void testM3DBlock();
    void testPLCCube();
    void testPLCCube2();
};
/*----------------------------------------------------------------------------*/
#endif /* TETGENFACADETESTSUITE_H_ */
/*----------------------------------------------------------------------------*/


/*----------------------------------------------------------------------------*/
/*
 * GeomToolKitTestSuite.h
 *
 *  Created on: 20 janv. 2011
 *      Author: ledouxf
 */
/*----------------------------------------------------------------------------*/
#ifndef GEOMTOOLKITTESTSUITE_H_
#define GEOMTOOLKITTESTSUITE_H_
/*----------------------------------------------------------------------------*/
#include <string>
#include <cppunit/extensions/HelperMacros.h>
/*----------------------------------------------------------------------------*/
class GeomToolKitTestSuite: public CppUnit::TestFixture {

    CPPUNIT_TEST_SUITE(GeomToolKitTestSuite);
    CPPUNIT_TEST(testInCircumCircle1);
    CPPUNIT_TEST(testInCircumCircle2);
    CPPUNIT_TEST(testInCircumCircle3);
    CPPUNIT_TEST(testInCircumCircle4);
    CPPUNIT_TEST(testInCircumCircle5);
    CPPUNIT_TEST(testInCircumCircle6);
    CPPUNIT_TEST(testInCircumCircle7);
    CPPUNIT_TEST(testInCircumCircle8);
    CPPUNIT_TEST(testInCircumCircle9);
    CPPUNIT_TEST(testInCircumCircle10);
    CPPUNIT_TEST(testInCircumCircle11);
    CPPUNIT_TEST(testInCircumCircle12);
    CPPUNIT_TEST_SUITE_END();

public:
    void setUp(){}

    void tearDown(){}

    void testInCircumCircle1();
    void testInCircumCircle2();
    void testInCircumCircle3();
    void testInCircumCircle4();
    void testInCircumCircle5();
    void testInCircumCircle6();
    void testInCircumCircle7();
    void testInCircumCircle8();
    void testInCircumCircle9();
    void testInCircumCircle10();
    void testInCircumCircle11();
    void testInCircumCircle12();

};
/*----------------------------------------------------------------------------*/
#endif /* GEOMTOOLKITTESTSUITE_H_ */
/*----------------------------------------------------------------------------*/




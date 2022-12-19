/*----------------------------------------------------------------------------*/
/** \file TritonUnitTestsLauncher
 *
 *  \author Franck Ledoux
 *
 *  \date 19/11/2011
 */
/*----------------------------------------------------------------------------*/
#include <cppunit/ui/text/TestRunner.h>
#include <cppunit/extensions/HelperMacros.h>
/*----------------------------------------------------------------------------*/
#include <iostream>
/*----------------------------------------------------------------------------*/
//#include "Delaunay2DTestSuite.h"
//#include "EdgeLength2DTestSuite.h"
//#include "GeomToolKitTestSuite.h"
//#include "RobustGeomToolKitTestSuite.h"
#include "TetgenFacadeTestSuite.h"
#include "NetgenFacadeTestSuite.h"
/*----------------------------------------------------------------------------*/
CPPUNIT_TEST_SUITE_REGISTRATION(NetgenFacadeTestSuite);
CPPUNIT_TEST_SUITE_REGISTRATION(TetgenFacadeTestSuite);
//CPPUNIT_TEST_SUITE_REGISTRATION(Delaunay2DTestSuite);
//CPPUNIT_TEST_SUITE_REGISTRATION(EdgeLength2DTestSuite);
//CPPUNIT_TEST_SUITE_REGISTRATION(GeomToolKitTestSuite);
//CPPUNIT_TEST_SUITE_REGISTRATION(RobustGeomToolKitTestSuite);

/*----------------------------------------------------------------------------*/
int main(int argc, char** argv){

    std::cout<<"##########################"<<"\n";
    std::cout<<"#### TRITON UNIT TEST ####"<<"\n";
    std::cout<<"##########################"<<"\n";
    CppUnit::TextUi::TestRunner runner;
    CppUnit::TestFactoryRegistry &registry =
                                CppUnit::TestFactoryRegistry::getRegistry();
    runner.addTest(registry.makeTest());
    runner.run();
    return 0;
}
/*----------------------------------------------------------------------------*/


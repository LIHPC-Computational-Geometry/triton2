cmake_minimum_required(VERSION 3.14.0)

project(TRITON2_SUITE)

#==============================================================================
#Versionnage du projet

set(TRITON2_SUITE_VERSION_MAJOR 1)
set(TRITON2_SUITE_VERSION_MINOR 0)
set(TRITON2_SUITE_VERSION_PATCH 1)

set(TRITON2_SUITE_VERSION
        ${TRITON2_SUITE_VERSION_MAJOR}.${TRITON2_SUITE_VERSION_MINOR}.${TRITON2_SUITE_VERSION_PATCH})


option(BUILD_SHARED_LIBS    ON)

add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/Triton2)


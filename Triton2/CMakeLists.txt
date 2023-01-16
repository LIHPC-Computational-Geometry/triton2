#==============================================================================
# GMDSCEA CORE MODULE
#==============================================================================
project(triton2 CXX)

#add_definitions(-DDLLEXPORT)

#==============================================================================
# CODE STRUCTURE
#==============================================================================

include_directories("${CMAKE_CURRENT_SOURCE_DIR}/inc")

set(TRITON2_SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src)

set(TRITON2_HEADER_FILES
	inc/Triton2/TetgenInterface/TetgenFacade.h
        inc/Triton2/TetgenInterface/tetgen.h
        )


set(TRITON2_SRC_FILES
        src/NetgenFacade.cpp
        src/TetgenFacade.cpp
	src/tetgen/tetgen.cxx
	src/tetgen/predicates.cxx
        )

set_source_files_properties(src/tetgen/predicates.cxx PROPERTIES COMPILE_FLAGS -O0)

#==============================================================================
# DEPENDENCIES
#==============================================================================
find_package(gmds REQUIRED)
include_directories(${GMDS_INCLUDE_DIR})

find_package(gmdscea REQUIRED)
include_directories(${GMDSCEA_INCLUDE_DIR})

#==============================================================================
# TEST MODULE
#==============================================================================

#set(GMDS_TEST ON CACHE BOOL "GMDS Kernel test module")
#if (GMDS_TEST)
#  add_subdirectory(test)
#endif ()

#==============================================================================
# LIBRARY BUILDING
#==============================================================================

add_library(triton2 ${TRITON2_HEADER_FILES} ${TRITON2_SRC_FILES})

#set_target_properties(gmds PROPERTIES PUBLIC_HEADER "${GMDS_HEADER_FILES}")


#target_compile_features(gmds PUBLIC cxx_std_11)


target_link_libraries(triton2 PRIVATE gmds gmdscea)

#add_custom_command(
#   TARGET gmds
#   POST_BUILD
#   COMMAND ${CMAKE_COMMAND} ARGS -E copy_directory  $<TARGET_FILE:gmds> ${CMAKE_BINARY_DIR}/gen/
#   COMMENT "Copy gmds lib"
#   VERBATIM
#   )

#target_include_directories(gmds INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>/generated>)
target_include_directories(triton2 PUBLIC "$<INSTALL_INTERFACE:include>")

install(TARGETS triton2
        EXPORT triton2 DESTINATION
        LIBRARY DESTINATION lib
        PUBLIC_HEADER DESTINATION include
        )

install(FILES
        ${TRITON2_HEADER_FILES} DESTINATION include/Triton2)


set(INCLUDE_INSTALL_DIR include/ )
set(LIB_INSTALL_DIR lib/ )
set(SYSCONFIG_INSTALL_DIR etc/gmds/ )

include(CMakePackageConfigHelpers)
configure_package_config_file(
        triton2Config.cmake.in ${CMAKE_CURRENT_BINARY_DIR}/triton2Config.cmake
        INSTALL_DESTINATION ${LIB_INSTALL_DIR}/triton2/cmake
        PATH_VARS INCLUDE_INSTALL_DIR SYSCONFIG_INSTALL_DIR
)

write_basic_package_version_file(
        ${CMAKE_CURRENT_BINARY_DIR}/triton2ConfigVersion.cmake
        VERSION ${TRITON2_SUITE_VERSION}
        COMPATIBILITY SameMajorVersion
)

install(FILES
        ${CMAKE_CURRENT_BINARY_DIR}/triton2Config.cmake
        ${CMAKE_CURRENT_BINARY_DIR}/triton2ConfigVersion.cmake
        DESTINATION ${LIB_INSTALL_DIR}/triton2/cmake
        )


install(EXPORT triton2 DESTINATION lib/triton2)

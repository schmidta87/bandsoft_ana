# CMakeLists for calibration
project(calibration)

# Set up
cmake_minimum_required (VERSION 3.0)
set(CMAKE_CXX_FLAGS "-std=c++11")
set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake_modules")
find_package(ROOT REQUIRED)


FILE(GLOB BAND_LIBS $ENV{BANDSOFT_TOOLS_DIR}/lib/libclashit_tree.* $ENV{BANDSOFT_TOOLS_DIR}/lib/libbandhit_tree.* $ENV{BANDSOFT_TOOLS_DIR}/lib/libtaghit_tree.* $ENV{BANDSOFT_TOOLS_DIR}/lib/libgenpart_tree.*)


set(INCLUDES ${ROOT_INCLUDE_DIR} ${CMAKE_SOURCE_DIR}/include $ENV{BANDSOFT_TOOLS_DIR}/include )

SET(CMAKE_SKIP_BUILD_RPATH  FALSE)
SET(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)
SET(CMAKE_INSTALL_RPATH "${CMAKE_SOURCE_DIR}/lib")
SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
LIST(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES "${CMAKE_SOURCE_DIR}/lib" isSystemDir)
IF("${isSystemDir}" STREQUAL "-1")
   SET(CMAKE_INSTALL_RPATH "${CMAKE_SOURCE_DIR}/lib")
ENDIF("${isSystemDir}" STREQUAL "-1")

# Add directories
#add_subdirectory(bandhit)
#add_subdirectory(clashit)
#add_subdirectory(taghit)
#add_subdirectory(genpart)
#add_subdirectory(doubleratio)
add_subdirectory(background)
add_subdirectory(inclusive)
add_subdirectory(efficiency)
add_subdirectory(final_skim)
add_subdirectory(photon_align)
add_subdirectory(bad_bars)
add_subdirectory(rcdb_info)

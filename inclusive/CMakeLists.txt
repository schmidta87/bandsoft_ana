set(INCLUDES ${ROOT_INCLUDE_DIR} ${CMAKE_SOURCE_DIR}/include )

include_directories(${INCLUDES} )

FILE(GLOB BAND_LIBS $ENV{BANDSOFT_TOOLS_DIR}/lib/libclas12fiducial.so)

set(CMAKE_VERBOSE_MAKEFILE ON)



add_executable(compare compare.cpp )
target_link_libraries(compare bandhit_tree clashit_tree taghit_tree ${ROOT_LIBRARIES} )
install( TARGETS compare
		DESTINATION ${CMAKE_SOURCE_DIR}/bin ) 

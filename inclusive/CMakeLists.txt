set(INCLUDES ${ROOT_INCLUDE_DIR} ${CMAKE_SOURCE_DIR}/include )

include_directories(${INCLUDES} )

FILE(GLOB BAND_LIBS $ENV{BANDSOFT_TOOLS_DIR}/lib/libclas12fiducial.so)

set(CMAKE_VERBOSE_MAKEFILE ON)



add_executable(compare compare.cpp )
target_link_libraries(compare bandhit_tree clashit_tree taghit_tree genpart_tree ${ROOT_LIBRARIES} )
install( TARGETS compare
		DESTINATION ${CMAKE_SOURCE_DIR}/bin ) 

add_executable(grid grid.cpp )
target_link_libraries(grid bandhit_tree clashit_tree taghit_tree genpart_tree ${ROOT_LIBRARIES} )
install( TARGETS grid
		DESTINATION ${CMAKE_SOURCE_DIR}/bin ) 

set(INCLUDES ${ROOT_INCLUDE_DIR} ${CMAKE_SOURCE_DIR}/include )

include_directories(${INCLUDES} )

FILE(GLOB BAND_LIBS $ENV{BANDSOFT_TOOLS_DIR}/lib/libclas12fiducial.so)

set(CMAKE_VERBOSE_MAKEFILE ON)


add_executable(grid_neutron grid_neutron.cpp )
target_link_libraries(grid_neutron bandhit_tree clashit_tree taghit_tree genpart_tree ${ROOT_LIBRARIES} )
install( TARGETS grid_neutron
		DESTINATION ${CMAKE_SOURCE_DIR}/bin ) 

add_executable(test_pointing test_pointing.cpp )
target_link_libraries(test_pointing bandhit_tree clashit_tree taghit_tree genpart_tree ${ROOT_LIBRARIES} )
install( TARGETS test_pointing
		DESTINATION ${CMAKE_SOURCE_DIR}/bin ) 

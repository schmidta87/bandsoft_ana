
include_directories(${INCLUDES} )


set(CMAKE_VERBOSE_MAKEFILE ON)

add_executable(grid_neutron grid_neutron.cpp )
target_link_libraries(grid_neutron ${ROOT_LIBRARIES} ${BAND_LIBS})
install( TARGETS grid_neutron
		DESTINATION ${CMAKE_SOURCE_DIR}/bin ) 

add_executable(test_pointing test_pointing.cpp )
target_link_libraries(test_pointing ${ROOT_LIBRARIES} ${BAND_LIBS})
install( TARGETS test_pointing
		DESTINATION ${CMAKE_SOURCE_DIR}/bin ) 

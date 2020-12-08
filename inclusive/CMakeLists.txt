
include_directories(${INCLUDES} )


set(CMAKE_VERBOSE_MAKEFILE ON)

add_executable(compare compare.cpp )
target_link_libraries(compare ${ROOT_LIBRARIES} ${BAND_LIBS})
install( TARGETS compare
		DESTINATION ${CMAKE_SOURCE_DIR}/bin ) 

add_executable(grid grid.cpp )
target_link_libraries(grid ${ROOT_LIBRARIES} ${BAND_LIBS})
install( TARGETS grid
		DESTINATION ${CMAKE_SOURCE_DIR}/bin ) 

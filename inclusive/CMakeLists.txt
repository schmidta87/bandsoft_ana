set(INCLUDES ${ROOT_INCLUDE_DIR} ${CMAKE_SOURCE_DIR}/include )

include_directories(${INCLUDES} )

FILE(GLOB BAND_LIBS $ENV{BANDSOFT_TOOLS_DIR}/lib/libclashit_tree.dylib $ENV{BANDSOFT_TOOLS_DIR}/lib/libbandhit_tree.dylib $ENV{BANDSOFT_TOOLS_DIR}/lib/libtaghit_tree.dylib $ENV{BANDSOFT_TOOLS_DIR}/lib/libgenpart_tree.dylib)

set(CMAKE_VERBOSE_MAKEFILE ON)

add_executable(compare compare.cpp )
target_link_libraries(compare ${ROOT_LIBRARIES} ${BAND_LIBS})
install( TARGETS compare
		DESTINATION ${CMAKE_SOURCE_DIR}/bin ) 

add_executable(grid grid.cpp )
target_link_libraries(grid ${ROOT_LIBRARIES} ${BAND_LIBS})
install( TARGETS grid
		DESTINATION ${CMAKE_SOURCE_DIR}/bin ) 

set(INCLUDES ${ROOT_INCLUDE_DIR} ${CMAKE_SOURCE_DIR}/include )

include_directories(${INCLUDES} )

FILE(GLOB BAND_LIBS $ENV{BANDSOFT_TOOLS_DIR}/lib/libclashit_tree.dylib $ENV{BANDSOFT_TOOLS_DIR}/lib/libbandhit_tree.dylib $ENV{BANDSOFT_TOOLS_DIR}/lib/libtaghit_tree.dylib $ENV{BANDSOFT_TOOLS_DIR}/lib/libgenpart_tree.dylib)

set(CMAKE_VERBOSE_MAKEFILE ON)

add_executable(final_inclusive final_inclusive.cpp )
target_link_libraries(final_inclusive ${ROOT_LIBRARIES} ${BAND_LIBS})
install( TARGETS final_inclusive
		DESTINATION ${CMAKE_SOURCE_DIR}/bin ) 

set(INCLUDES ${ROOT_INCLUDE_DIR} ${CMAKE_SOURCE_DIR}/include )

include_directories(${INCLUDES} )

FILE(GLOB BAND_LIBS $ENV{BANDSOFT_TOOLS_DIR}/lib/libclashit_tree.dylib $ENV{BANDSOFT_TOOLS_DIR}/lib/libbandhit_tree.dylib $ENV{BANDSOFT_TOOLS_DIR}/lib/libtaghit_tree.dylib $ENV{BANDSOFT_TOOLS_DIR}/lib/libgenpart_tree.dylib)

set(CMAKE_VERBOSE_MAKEFILE ON)

add_executable(test_skim test_skim.cpp )
target_link_libraries(test_skim ${ROOT_LIBRARIES} ${BAND_LIBS})
install( TARGETS test_skim
		DESTINATION ${CMAKE_SOURCE_DIR}/bin ) 

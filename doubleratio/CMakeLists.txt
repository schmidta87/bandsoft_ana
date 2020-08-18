set(INCLUDES ${ROOT_INCLUDE_DIR} ${CMAKE_SOURCE_DIR}/include )
include_directories(${INCLUDES})

set(CMAKE_VERBOSE_MAKEFILE ON)

add_executable(v0 v0.cpp )
target_link_libraries(v0 bandhit_tree clashit_tree taghit_tree ${ROOT_LIBRARIES} )
install( TARGETS v0
		DESTINATION ${CMAKE_SOURCE_DIR}/bin )

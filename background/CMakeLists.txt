set(INCLUDES ${ROOT_INCLUDE_DIR} ${CMAKE_SOURCE_DIR}/include )
include_directories(${INCLUDES})

set(CMAKE_VERBOSE_MAKEFILE ON)

add_executable(event_mixing event_mixing.cpp )
target_link_libraries(event_mixing bandhit_tree clashit_tree taghit_tree ${ROOT_LIBRARIES} )
install( TARGETS event_mixing
		DESTINATION ${CMAKE_SOURCE_DIR}/bin ) 


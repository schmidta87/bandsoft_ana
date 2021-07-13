
include_directories(${INCLUDES})


set(CMAKE_VERBOSE_MAKEFILE ON)

add_executable(compare_xdist compare_xdist.cpp )
target_link_libraries(compare_xdist ${ROOT_LIBRARIES} ${BAND_LIBS})
install( TARGETS compare_xdist
		DESTINATION ${CMAKE_SOURCE_DIR}/bin )

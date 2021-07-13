
include_directories(${INCLUDES} )


set(CMAKE_VERBOSE_MAKEFILE ON)

add_executable(get_gated get_gated.cpp )
target_link_libraries(get_gated ${ROOT_LIBRARIES} ${BAND_LIBS})
install( TARGETS get_gated
		DESTINATION ${CMAKE_SOURCE_DIR}/bin ) 


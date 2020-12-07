set(INCLUDES ${ROOT_INCLUDE_DIR} ${CMAKE_SOURCE_DIR}/include )

include_directories(${INCLUDES} )


set(CMAKE_VERBOSE_MAKEFILE ON)

add_executable(final_inclusive final_inclusive.cpp )
target_link_libraries(final_inclusive ${ROOT_LIBRARIES} ${BAND_LIBS})
install( TARGETS final_inclusive
		DESTINATION ${CMAKE_SOURCE_DIR}/bin ) 

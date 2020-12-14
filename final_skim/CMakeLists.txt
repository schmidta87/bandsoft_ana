
include_directories(${INCLUDES} )


set(CMAKE_VERBOSE_MAKEFILE ON)

add_executable(final_inclusive final_inclusive.cpp )
target_link_libraries(final_inclusive ${ROOT_LIBRARIES} ${BAND_LIBS})
install( TARGETS final_inclusive
		DESTINATION ${CMAKE_SOURCE_DIR}/bin ) 

add_executable(final_tagged final_tagged.cpp )
target_link_libraries(final_tagged ${ROOT_LIBRARIES} ${BAND_LIBS})
install( TARGETS final_tagged
		DESTINATION ${CMAKE_SOURCE_DIR}/bin ) 

add_executable(final_neutrons final_neutrons.cpp )
target_link_libraries(final_neutrons ${ROOT_LIBRARIES} ${BAND_LIBS})
install( TARGETS final_neutrons
		DESTINATION ${CMAKE_SOURCE_DIR}/bin ) 


include_directories(${INCLUDES} )


set(CMAKE_VERBOSE_MAKEFILE ON)

add_executable(final_inclusive final_inclusive.cpp )
target_link_libraries(final_inclusive ${ROOT_LIBRARIES} ${BAND_LIBS})
install( TARGETS final_inclusive
		DESTINATION ${CMAKE_SOURCE_DIR}/bin ) 

add_executable(final_tagged_newnorm final_tagged_newnorm.cpp )
target_link_libraries(final_tagged_newnorm ${ROOT_LIBRARIES} ${BAND_LIBS})
install( TARGETS final_tagged_newnorm
		DESTINATION ${CMAKE_SOURCE_DIR}/bin ) 

add_executable(final_neutrons final_neutrons.cpp )
target_link_libraries(final_neutrons ${ROOT_LIBRARIES} ${BAND_LIBS})
install( TARGETS final_neutrons
		DESTINATION ${CMAKE_SOURCE_DIR}/bin ) 

add_executable(get_normalization get_normalization.cpp )
target_link_libraries(get_normalization ${ROOT_LIBRARIES} ${BAND_LIBS})
install( TARGETS get_normalization
		DESTINATION ${CMAKE_SOURCE_DIR}/bin ) 

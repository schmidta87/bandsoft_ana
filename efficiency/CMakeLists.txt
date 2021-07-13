include_directories(${INCLUDES} )


add_executable(calc_eff calc_eff.cpp )
target_link_libraries(calc_eff ${ROOT_LIBRARIES} ${BAND_LIBS})
install( TARGETS calc_eff
		DESTINATION ${CMAKE_SOURCE_DIR}/bin ) 


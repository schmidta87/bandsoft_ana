include_directories(${INCLUDES} )


add_executable(resolution resolution.cpp )
target_link_libraries(resolution ${ROOT_LIBRARIES} ${BAND_LIBS})
install( TARGETS resolution
		DESTINATION ${CMAKE_SOURCE_DIR}/bin/data_sim ) 

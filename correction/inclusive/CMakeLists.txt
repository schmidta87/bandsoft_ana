
include_directories(${INCLUDES} )

add_executable(raw_yield_inclusive raw_yield_inclusive.cpp )
target_link_libraries(raw_yield_inclusive ${ROOT_LIBRARIES} ${BAND_LIBS})
install( TARGETS raw_yield_inclusive
		DESTINATION ${CMAKE_SOURCE_DIR}/bin/cross_section ) 


add_executable(bin_migration_inclusive bin_migration_inclusive.cpp )
target_link_libraries(bin_migration_inclusive ${ROOT_LIBRARIES} ${BAND_LIBS})
install( TARGETS bin_migration_inclusive
		DESTINATION ${CMAKE_SOURCE_DIR}/bin/cross_section ) 

	#add_executable(radiative_inclusive radiative_inclusive.cpp )
	#target_link_libraries(radiative_inclusive ${ROOT_LIBRARIES} ${BAND_LIBS})
	#install( TARGETS cross_section/radiative_inclusive
	#		DESTINATION ${CMAKE_SOURCE_DIR}/bin ) 

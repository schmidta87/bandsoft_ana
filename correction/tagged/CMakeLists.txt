
include_directories(${INCLUDES} )

add_executable(raw_yield_tagged raw_yield_tagged.cpp )
target_link_libraries(raw_yield_tagged ${ROOT_LIBRARIES} ${BAND_LIBS})
install( TARGETS raw_yield_tagged
		DESTINATION ${CMAKE_SOURCE_DIR}/bin/cross_section ) 


add_executable(bin_migration_tagged bin_migration_tagged.cpp )
target_link_libraries(bin_migration_tagged ${ROOT_LIBRARIES} ${BAND_LIBS})
install( TARGETS bin_migration_tagged
		DESTINATION ${CMAKE_SOURCE_DIR}/bin/cross_section ) 

add_executable(radiative_tagged radiative_tagged.cpp )
target_link_libraries(radiative_tagged ${ROOT_LIBRARIES} ${BAND_LIBS})
install( TARGETS radiative_tagged
		DESTINATION ${CMAKE_SOURCE_DIR}/bin/cross_section ) 

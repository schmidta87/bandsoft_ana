
include_directories(${INCLUDES} )

add_executable(raw_yield raw_yield.cpp )
target_link_libraries(raw_yield ${ROOT_LIBRARIES} ${BAND_LIBS})
install( TARGETS raw_yield
		DESTINATION ${CMAKE_SOURCE_DIR}/bin ) 


add_executable(bin_migration bin_migration.cpp )
target_link_libraries(bin_migration ${ROOT_LIBRARIES} ${BAND_LIBS})
install( TARGETS bin_migration
		DESTINATION ${CMAKE_SOURCE_DIR}/bin ) 

add_executable(radiative radiative.cpp )
target_link_libraries(radiative ${ROOT_LIBRARIES} ${BAND_LIBS})
install( TARGETS radiative
		DESTINATION ${CMAKE_SOURCE_DIR}/bin ) 

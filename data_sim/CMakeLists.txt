
include_directories(${INCLUDES} )

add_executable(yield_ratio yield_ratio.cpp )
target_link_libraries(yield_ratio ${ROOT_LIBRARIES} ${BAND_LIBS})
install( TARGETS yield_ratio
		DESTINATION ${CMAKE_SOURCE_DIR}/bin/data_sim ) 
add_executable(yield_ratio_virt yield_ratio_virt.cpp )
target_link_libraries(yield_ratio_virt ${ROOT_LIBRARIES} ${BAND_LIBS})
install( TARGETS yield_ratio_virt
		DESTINATION ${CMAKE_SOURCE_DIR}/bin/data_sim ) 


include_directories(${INCLUDES} )


add_executable(nalgo_multiplicities nalgo_multiplicities.cpp )
target_link_libraries(nalgo_multiplicities ${ROOT_LIBRARIES} ${BAND_LIBS})
install( TARGETS nalgo_multiplicities
		DESTINATION ${CMAKE_SOURCE_DIR}/bin ) 


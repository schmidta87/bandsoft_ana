
include_directories(${INCLUDES})


set(CMAKE_VERBOSE_MAKEFILE ON)

#add_executable(event_mixing event_mixing.cpp )
#target_link_libraries(event_mixing ${ROOT_LIBRARIES} ${BAND_LIBS})
#install( TARGETS event_mixing
#		DESTINATION ${CMAKE_SOURCE_DIR}/bin ) 

add_executable(event_mixing_tof event_mixing_tof.cpp )
target_link_libraries(event_mixing_tof ${ROOT_LIBRARIES} ${BAND_LIBS})
install( TARGETS event_mixing_tof
		DESTINATION ${CMAKE_SOURCE_DIR}/bin ) 


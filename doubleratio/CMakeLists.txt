
include_directories(${INCLUDES} )


set(CMAKE_VERBOSE_MAKEFILE ON)

add_executable(v0 v0.cpp )
target_link_libraries(v0 ${ROOT_LIBRARIES} ${BAND_LIBS})
install( TARGETS v0
		DESTINATION ${CMAKE_SOURCE_DIR}/bin ) 

add_executable(v0_noBG v0_noBG.cpp )
target_link_libraries(v0_noBG ${ROOT_LIBRARIES} ${BAND_LIBS})
install( TARGETS v0_noBG
		DESTINATION ${CMAKE_SOURCE_DIR}/bin )
	
	#add_executable(v0_noBG_noTClone v0_noBG_noTClone.cpp )
	#target_link_libraries(v0_noBG ${ROOT_LIBRARIES} ${BAND_LIBS})
	#install( TARGETS v0_noBG_noTClone
	#DESTINATION ${CMAKE_SOURCE_DIR}/bin )

add_executable(good_events good_events.cpp )
target_link_libraries(good_events ${ROOT_LIBRARIES} ${BAND_LIBS})
install( TARGETS good_events
		DESTINATION ${CMAKE_SOURCE_DIR}/bin )

add_executable(good_events_noBG good_events_noBG.cpp )
target_link_libraries(good_events_noBG ${ROOT_LIBRARIES} ${BAND_LIBS})
install( TARGETS good_events_noBG
		DESTINATION ${CMAKE_SOURCE_DIR}/bin )

add_executable(good_events_noBG_noTClone good_events_noBG_noTClone.cpp )
target_link_libraries(good_events_noBG_noTClone ${ROOT_LIBRARIES} ${BAND_LIBS})
install( TARGETS good_events_noBG_noTClone
		DESTINATION ${CMAKE_SOURCE_DIR}/bin )

add_executable(make_histograms make_histograms.cpp )
target_link_libraries(make_histograms ${ROOT_LIBRARIES} ${BAND_LIBS})
install( TARGETS make_histograms
		DESTINATION ${CMAKE_SOURCE_DIR}/bin )





include_directories(${INCLUDES})


set(CMAKE_VERBOSE_MAKEFILE ON)

add_executable(initial_bar initial_bar.cpp )
target_link_libraries(initial_bar ${ROOT_LIBRARIES} ${BAND_LIBS})
install( TARGETS initial_bar
		DESTINATION ${CMAKE_SOURCE_DIR}/bin )

add_executable(initial_run initial_run.cpp calib_helper.cpp )
target_link_libraries(initial_run ${ROOT_LIBRARIES} ${BAND_LIBS})
install( TARGETS initial_run
		DESTINATION ${CMAKE_SOURCE_DIR}/bin )

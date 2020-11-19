set(INCLUDES ${ROOT_INCLUDE_DIR} ${CMAKE_SOURCE_DIR}/include )

include_directories(${INCLUDES} $ENV{BANDSOFT_TOOLS_DIR}/include)

FILE(GLOB BAND_LIBS $ENV{BANDSOFT_TOOLS_DIR}/lib/libclas12fiducial.so)

set(CMAKE_VERBOSE_MAKEFILE ON)

add_executable(v0 v0.cpp )
target_link_libraries(v0 bandhit_tree clashit_tree taghit_tree ${ROOT_LIBRARIES} ${BAND_LIBS})
install( TARGETS v0
		DESTINATION ${CMAKE_SOURCE_DIR}/bin ) 

add_executable(v0_noBG v0_noBG.cpp )
target_link_libraries(v0_noBG bandhit_tree clashit_tree taghit_tree ${ROOT_LIBRARIES} ${BAND_LIBS})
install( TARGETS v0_noBG
		DESTINATION ${CMAKE_SOURCE_DIR}/bin )

add_executable(v0_noBG_noTClone v0_noBG_noTClone.cpp )
target_link_libraries(v0_noBG_noTClone bandhit_tree clashit_tree taghit_tree ${ROOT_LIBRARIES} ${BAND_LIBS})
install( TARGETS v0_noBG_noTClone
		DESTINATION ${CMAKE_SOURCE_DIR}/bin )

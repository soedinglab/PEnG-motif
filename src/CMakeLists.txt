include_directories(shared)
add_subdirectory(shared)

include_directories("${PROJECT_BINARY_DIR}")

set(PENG_SOURCE
	${PENG_SOURCE}
	${CMAKE_CURRENT_SOURCE_DIR}/Global.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/Global.h
	${CMAKE_CURRENT_SOURCE_DIR}/PatternCensus.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/PatternCensus.h
	${CMAKE_CURRENT_SOURCE_DIR}/log.h
)

add_library(PENG_OBJECTS OBJECT ${PENG_SOURCE})

# add the executable
add_executable(peng_motif 
	main.cpp	$<TARGET_OBJECTS:PENG_OBJECTS>
)
target_link_libraries (peng_motif LINK_PUBLIC shared)

INSTALL(PROGRAMS 
  ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/peng_motif
  DESTINATION bin
)

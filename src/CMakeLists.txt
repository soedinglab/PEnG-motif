include_directories(shared)
add_subdirectory(shared)

include_directories("${PROJECT_BINARY_DIR}")

set(PENG_SOURCE
	${PENG_SOURCE}
	${CMAKE_CURRENT_SOURCE_DIR}/Global.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/Global.h
	${CMAKE_CURRENT_SOURCE_DIR}/peng.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/peng.h
	${CMAKE_CURRENT_SOURCE_DIR}/log.h
	${CMAKE_CURRENT_SOURCE_DIR}/iupac_alphabet.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/iupac_alphabet.h
	${CMAKE_CURRENT_SOURCE_DIR}/iupac_pattern.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/iupac_pattern.h
	${CMAKE_CURRENT_SOURCE_DIR}/base_pattern.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/base_pattern.h
	${CMAKE_CURRENT_SOURCE_DIR}/helper-inl.h
		utils.h)

add_library(PENG_OBJECTS OBJECT ${PENG_SOURCE})

# add the executable
add_executable(peng_motif 
	main.cpp	$<TARGET_OBJECTS:PENG_OBJECTS>)
target_link_libraries (peng_motif LINK_PUBLIC shared)

INSTALL(PROGRAMS 
  ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/peng_motif
  DESTINATION bin
)

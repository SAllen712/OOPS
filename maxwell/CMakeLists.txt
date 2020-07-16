cmake_minimum_required(VERSION 3.0)
project(Maxwell)

#Generate Parameter and Parser files
set(PARAM_SRC ${CMAKE_CURRENT_SOURCE_DIR}/src/maxwellparser.cpp)
set(PARAM_INC
	${CMAKE_CURRENT_SOURCE_DIR}/include/maxwellparameters.h
	${CMAKE_CURRENT_SOURCE_DIR}/include/maxwellparser.h
	)
set(SETUP_SRC ${CMAKE_CURRENT_SOURCE_DIR}/scripts/maxwell.json)

add_custom_command(
	OUTPUT ${PARAM_INC}
	       ${PARAM_SRC}
	DEPENDS ${SETUP_SRC}
	COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_SOURCE_DIR}/scripts/genParams.py ${SETUP_SRC}
	COMMENT "Generating custom Parameters files"
	WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
	VERBATIM USES_TERMINAL
)

set(MAXWELL_INCLUDE_FILES
	include/maxwell.h
	)
set(MAXWELL_SOURCE_FILES
	src/maxwell.cpp
	src/main.cpp
	)

set(SOURCE_FILES ${MAXWELL_INCLUDE_FILES} ${MAXWELL_SOURCE_FILES} ${PARAM_INC} ${PARAM_SRC})
add_executable(Maxwell ${SOURCE_FILES})
target_include_directories(Maxwell PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_include_directories(Maxwell PRIVATE ${CMAKE_SOURCE_DIR}/include)
target_link_libraries(Maxwell oops ${EXTRA_LIBS})


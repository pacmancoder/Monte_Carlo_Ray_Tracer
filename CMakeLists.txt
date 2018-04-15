cmake_minimum_required(VERSION 3.8)

project(mcrt)

# Set gcc / g++ compiler to be able to use openMP 
SET(CMAKE_C_COMPILER "gcc")
SET(CMAKE_CXX_COMPILER "g++")

FIND_PACKAGE( OpenMP )
if(OPENMP_FOUND)
	message("OPENMP FOUND")
	set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
	set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

include_directories(${OPENMP_INCLUDE_DIR})

# Our own include directory
set(COMMON_INCLUDES ${PROJECT_SOURCE_DIR}/include)
include_directories(${COMMON_INCLUDES})

# Source files
file(GLOB INTERNAL_SOURCE
        ${PROJECT_SOURCE_DIR}/src/*.cpp
        ${PROJECT_SOURCE_DIR}/external_libraries/common_src/*.cpp)
file(GLOB INTERNAL_HEADERS
        ${PROJECT_SOURCE_DIR}/include/*.h
        ${PROJECT_SOURCE_DIR}/external_libraries/common_include/*.h)

# Create the executable from our sources
add_executable(${PROJECT_NAME} ${INTERNAL_SOURCE} ${INTERNAL_HEADERS})

target_link_libraries(${PROJECT_NAME} PUBLIC glm)
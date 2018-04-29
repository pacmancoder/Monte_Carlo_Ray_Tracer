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


include_directories(libs/include)
include_directories(include)

# Source files
file(GLOB INTERNAL_SOURCE
		libs/src/*.cpp
		src/*.cpp)

file(GLOB INTERNAL_HEADERS
        include/mcrt/*.h
        libs/include/external/*.h
        libs/include/external/kdtree++/*.hpp)

# Create the executable from our sources
add_library(${PROJECT_NAME} STATIC ${INTERNAL_SOURCE} ${INTERNAL_HEADERS})

target_link_libraries(${PROJECT_NAME} PUBLIC glm)

target_include_directories(${PROJECT_NAME} PUBLIC libs/include)
target_include_directories(${PROJECT_NAME} PUBLIC include)
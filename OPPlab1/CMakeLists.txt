cmake_minimum_required(VERSION 3.14)

set(CMAKE_CXX_STANDARD 17)

project(OPPlab1)

set(exe_name ${PROJECT_NAME})

set(MPI_ROOT "C:/Program Files (x86)/Microsoft SDKs/MPI")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/out)

include_directories(${MPI_ROOT}/include)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)
link_directories(${MPI_ROOT}/Lib/x64)

add_executable(
        ${exe_name} src/main.cpp
        src/Matrix.h
        src/Matrix.cpp
        src/Solution.h
        src/UsualSolution.cpp
        src/FirstSolution.cpp
        src/SecondSolution.cpp
)
target_link_libraries(${exe_name} msmpi.lib)

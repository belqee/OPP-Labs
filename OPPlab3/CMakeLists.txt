cmake_minimum_required(VERSION 3.14)

set(CMAKE_CXX_STANDARD 17)

project(OPPlab3)

set(exe_name ${PROJECT_NAME})

set(MPI_ROOT "C:/Program Files (x86)/Microsoft SDKs/MPI")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/out)

include_directories(${MPI_ROOT}/include)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)
link_directories(${MPI_ROOT}/Lib/x64)

add_executable(
        ${exe_name} src/main.c

)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /O3 /Wall /WX")

target_link_libraries(${exe_name} msmpi.lib)

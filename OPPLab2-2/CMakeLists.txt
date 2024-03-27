cmake_minimum_required(VERSION 3.14)

set(CMAKE_CXX_STANDARD 17)

project(OPPlab2)

set(exe_name ${PROJECT_NAME})

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/out)


include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)


add_executable(
        ${exe_name} src/main.cpp
        src/Matrix.h


)
add_compile_options(/O3 /Wall /WX)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
target_link_options(${exe_name} PRIVATE "-fopenmp")
# Добавляем флаги компиляции для OpenMP
find_package(OpenMP REQUIRED)
if(OpenMP_CXX_FOUND)
    target_link_libraries(${exe_name} PUBLIC OpenMP::OpenMP_CXX)
endif()


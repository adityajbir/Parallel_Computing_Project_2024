cmake_minimum_required(VERSION 3.12)
find_package(MPI REQUIRED)
find_package(caliper REQUIRED)
find_package(adiak REQUIRED)

message(STATUS "MPI includes : ${MPI_INCLUDE_PATH}")
message(STATUS "Caliper includes : ${caliper_INCLUDE_DIR}")
message(STATUS "Adiak includes : ${adiak_INCLUDE_DIRS}")

include_directories(SYSTEM ${MPI_INCLUDE_PATH})
include_directories(${caliper_INCLUDE_DIR})
include_directories(${adiak_INCLUDE_DIRS})
include_directories(${CMAKE_SOURCE_DIR}/common)
include_directories(${CMAKE_SOURCE_DIR}/src)

# Collect all .cpp files from src and common directories
file(GLOB SRC_FILES ${CMAKE_SOURCE_DIR}/../src/*.cpp)
file(GLOB COMMON_FILES ${CMAKE_SOURCE_DIR}/../common/*.cpp)

add_executable(main_executable ${SRC_FILES} ${COMMON_FILES})
target_link_libraries(main_executable PRIVATE MPI::MPI_CXX)
target_link_libraries(main_executable PRIVATE caliper)
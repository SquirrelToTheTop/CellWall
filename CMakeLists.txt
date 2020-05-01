cmake_minimum_required(VERSION 3.10)

# project name
project(CellWall VERSION 2.0)

set(cellWall_VERSION_MAJOR 2)
set(cellWall_VERSION_MINOR 0)

# include header directory
set(CellWall_h_dir ${PROJECT_SOURCE_DIR}/headers)
include_directories(${CellWall_h_dir})

#include_directories("./headers")

message("Project bin dir : " ${PROJECT_BINARY_DIR})
message("Project source dir : " ${PROJECT_SOURCE_DIR})
message("Project header dir : " ${CellWall_h_dir})

find_package(OpenMP)

# config file
configure_file(${CellWall_h_dir}/cellWallConfig.h.in ${CellWall_h_dir}/cellWallConfig.h)

# add code
add_library(${PROJECT_NAME}_obj SHARED ${PROJECT_SOURCE_DIR}/src/cellWallUtils.cpp
                                       ${PROJECT_SOURCE_DIR}/src/cellWallIO.cpp
				       ${PROJECT_SOURCE_DIR}/src/cellWallForces.cpp
)

target_compile_options(${PROJECT_NAME}_obj PRIVATE -O3 -Wall -g -fPIC ${OpenMP_C_FLAGS} -lm)
target_link_libraries(${PROJECT_NAME}_obj -lpng -lz ${OpenMP_C_LIB_NAMES})

# build exec from main.cpp
add_executable(${PROJECT_NAME} src/main.cpp)
target_compile_options(${PROJECT_NAME} PRIVATE -O3 -Wall -g ${OpenMP_C_FLAGS})
target_link_libraries(${PROJECT_NAME} ${PROJECT_NAME}_obj -lpng -lz -lm ${OpenMP_C_LIB_NAMES})


cmake_minimum_required(VERSION 3.10)

# project name
project(CellWall VERSION 2.0)

set(cellWall_VERSION_MAJOR 2)
set(cellWall_VERSION_MINOR 0)

option(BUILD_DEBUG "Build in debug mode" OFF)

if(BUILD_DEBUG)
  add_definitions(-DDEBUG)
endif(BUILD_DEBUG)

# include header directory
set(CellWall_h_dir ${PROJECT_SOURCE_DIR}/headers)
include_directories(${CellWall_h_dir})

#include_directories("./headers")

message("Project bin dir : " ${PROJECT_BINARY_DIR})
message("Project source dir : " ${PROJECT_SOURCE_DIR})
message("Project header dir : " ${CellWall_h_dir})

find_package(MPI)

if(${MPI_FOUND})
  # include_directories()
  message("MPI headers: " ${MPI_CXX_INCLUDE_DIRS})
  message("MPI Link Flags: " ${MPI_CXX_LINK_FLAGS})
  message("MPI Libraries: " ${MPI_CXX_LIBRARIES})
  message("MPI Compile Op: " ${MPI_CXX_COMPILE_OPTIONS})
  message("MPI dir: " ${MPI_INCLUDE_PATH})

  if( NOT EXISTS ${MPI_CXX_INCUDE_DIRS})
    message("> ERROR Header directory does not exist !")

    # tres very mauvais mais bon si MPI faisait son taff
    include_directories(/usr/lib/x86_64-linux-gnu/openmpi/include/)
    message("> Manual Correction : /usr/lib/x86_64-linux-gnu/openmpi/include/")
  endif(NOT EXISTS ${MPI_CXX_INCUDE_DIRS})

endif(${MPI_FOUND})

find_package(OpenMP)

message("OpenMP flags : " ${OpenMP_CXX_FLAGS})

include_directories(${MPI_CXX_INCLUDE_DIRS})

# config file
configure_file(${CellWall_h_dir}/cellWallConfig.h.in ${CellWall_h_dir}/cellWallConfig.h)

# add code
add_library(${PROJECT_NAME}_obj SHARED ${PROJECT_SOURCE_DIR}/src/cellWallUtils.cpp
                                       ${PROJECT_SOURCE_DIR}/src/cellWallIO.cpp
                                       ${PROJECT_SOURCE_DIR}/src/cellWallForces.cpp
                                       ${PROJECT_SOURCE_DIR}/src/cellWallMonolayer.cpp
                                       ${PROJECT_SOURCE_DIR}/src/cellWallLipidLayer.cpp
                                       ${PROJECT_SOURCE_DIR}/src/cellWallOptimization.cpp
                                       ${PROJECT_SOURCE_DIR}/src/cellWallDebug.cpp
)

target_compile_options(${PROJECT_NAME}_obj PRIVATE -O3 
                                                   -Wall
                                                   -g
                                                   -fPIC
                                                   -lm 
                                                   ${OpenMP_CXX_FLAGS}
                                                   ${MPI_CXX_COMPILE_OPTIONS}
)
target_link_libraries(${PROJECT_NAME}_obj -lpng
                                          -lz 
                                          ${OpenMP_CXX_LIB_NAMES}
                                          ${MPI_CXX_LIBRARIES}                                          
)

# build exec from main.cpp
add_executable(${PROJECT_NAME} src/main.cpp)
target_compile_options(${PROJECT_NAME} PRIVATE -O3 -Wall -g 
                                               ${OpenMP_CXX_FLAGS}
                                               ${MPI_CXX_COMPILE_OPTIONS}
)

target_link_libraries(${PROJECT_NAME} ${PROJECT_NAME}_obj -lpng -lz -lm 
                                                           ${OpenMP_CXX_LIB_NAMES}
                                                           ${MPI_CXX_LIBRARIES}
)


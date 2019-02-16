cmake_minimum_required (VERSION 3.5)

# Executable name: 
##############################################################
SET(EXEC OUTFile)
##############################################################

project(${EXEC})

#include
include_directories(/usr/local/include/eigen3)
include_directories(/usr/local/include)
include_directories(/home/aniruddha/Desktop/Composite_Layup/CUDA_experiments/include)
include_directories(/home/aniruddha/Desktop/Composite_Layup/CUDA_experiments/data)
include_directories(${CUDA_INCLUDE_DIRS})

find_package(CUDA REQUIRED)

# support for Libnabo  
SET(CMAKE_CXX_FLAGS "-std=c++11 -O3")

set(nlopt_ADDED_LIBS /usr/local/lib/libnlopt.so)

set(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS} -gencode arch=compute_20,code=sm_20; 
	--std=c++11; CUDA_ARCH_BIN=3.5;
	BUILD_SHARED_LIBS=off; 
	CMAKE_CXX_FLAGS=-fPIC; 
	WITH_QT=off; 
	WITH_VTK=off; 
	WITH_GTK=off; 
	WITH_OPENGL=off)
list(APPEND CUDA_NVCC_FLAGS "--expt-relaxed-constexpr")
set(CUDA_NVCC_FLAGS -I/opt/intel/mkl/include)



set(UTLIB src/utilities.cpp 
	src/file_rw.cpp 
	src/transformation_utilities.cpp
	src/opt_obj_eigen.cpp)
##############################################################
set(MAIN_FILE src/test_file_eigen.cu)
##############################################################
cuda_add_executable(${EXEC} ${MAIN_FILE} ${UTLIB} src/opt_obj_eigen.cpp)


# set(UTLIB src/utilities.cpp 
# 	src/file_rw.cpp 
# 	src/transformation_utilities.cpp
# 	src/opt_obj.cu)
# ##############################################################
# set(MAIN_FILE src/test_file_cuda.cu)
# ##############################################################
# cuda_add_executable(${EXEC} ${MAIN_FILE} ${UTLIB} src/opt_obj.cu)


##LINKED LIBRARIES
#link library for nlopt 
target_link_libraries (${EXEC} ${nlopt_ADDED_LIBS})

set_target_properties(${EXEC} PROPERTIES LINKER_LANGUAGE CXX)
set_property(TARGET ${EXEC} PROPERTY CXX_STANDARD 14)

#Last Update: 21.08.2016
#author: Julien Wulf
#University: UZH
#E-Mail: jwulf@physik.uzh.ch
#Description:
#The macros to find the packages (Root, etc) are in the "CMake" folder.
#Every Folder contains a CMakeLists.txt, which defines the linking etc
#For Every Subsystem a Libary will be created 


#General Stuff (Paths)
cmake_minimum_required(VERSION 2.8)
set (CMAKE_AUTOMOC OFF)
set (CMAKE_C_COMPILER mpicc)
set (CMAKE_CXX_COMPILER mpicxx)
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

project(XSLMY_ANALYSIS)
set (CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/CMake" ${CMAKE_MODULE_PATH})
set (CMAKE_BINARY_DIR ${CMAKE_SOURCE_DIR}/build)
set (EXECUTABLE_OUTPUT_PATH ${CMAKE_SOURCE_DIR}/bin)
set (LIBRARY_OUTPUT_PATH ${CMAKE_BINARY_DIR}/lib)
set (CMAKE_VERBOSE_MAKEFILE NO) 

#All headers of the subsystems muss be placed here 
set (SOURCE ${CMAKE_CURRENT_SOURCE_DIR})
set (MY_ANALYSIS_INCLUDE_DIRS
	${SOURCE}/include/
)

#Some options for ccmake 
OPTION (MY_ANALYSIS_USE_MPI "Acclerate  the calculations with Multithreading" OFF)
OPTION (MY_ANALYSIS_VERBOSE "Verbose of the compiler output" OFF)

if(MY_ANALYSIS_VERBOSE)
	set(CMAKE_VERBOSE_MAKEFILE YES) 
endif(MY_ANALYSIS_VERBOSE)

#Find external libaries and headers
find_package(ROOT REQUIRED)

if (XSLMY_ANALYSIS_USE_MPI)
	find_package(MPI)
	set (CMAKE_C_COMPILER mpicc)
	set (CMAKE_CXX_COMPILER mpicxx)
	set (CMAKE_EXECUTABLE_SUFFIX _mpi) 	
	add_definitions (-DXSLMY_ANALYSIS_MPI)
endif(XSLMY_ANALYSIS_USE_MPI)


#Global Headers
include_directories(${MPI_INCLUDE_PATH} ${MY_ANALYSIS_INCLUDE_DIRS} ${ROOT_INCLUDE_DIRS})
add_subdirectory(MyAnalysis)

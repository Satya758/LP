cmake_minimum_required(VERSION 3.1)

project(linearprogramming)

set(CMAKE_BUILD_TYPE "Debug")
# set(CMAKE_BUILD_TYPE "Release")
set(CMAKE_CXX_FLAGS "-std=c++14 -fopenmp")

add_subdirectory(include)
add_subdirectory(src)

####### Required Libraries ########################
# Header libs
set(EIGEN_INC "" CACHE PATH "Eigen Include Directory") #No default value
set(BOOST_INC "" CACHE PATH "Boost Lib Include Directory") #No default value
set(TBB_INC "" CACHE PATH "TBB Lib Include Directory") #No default value
set(SS_INC "" CACHE PATH "SuiteSparse Lib Include Directory") #No default value
set(PASTIX_INC "" CACHE PATH "Pastix Lib Include Directory") #No default value

# Header libs
# Binaries
set(BOOST_LOG_LIB "" CACHE PATH "Boost log library") #No default value
set(SS_CHOLMOD_LIB "" CACHE PATH "SuiteSparse CHOLMOD") #No default value
set(METIS_LIB "" CACHE PATH "Metis location")
set(LAPACK_LIB "" CACHE PATH "Lapack location")
set(BLAS_LIB "" CACHE PATH "Blas location")
set(TBB_LIB "" CACHE PATH "TBB location")
set(PASTIX_LIB "" CACHE PATH "Pastix location")
set(SCOTCH_LIB "" CACHE PATH "Pastix location")
# Binaries
####### Required Libraries ########################

include_directories(${EIGEN_INC})
include_directories(${BOOST_INC})
include_directories(${SS_INC})
include_directories(${TBB_INC})
include_directories(${PASTIX_INC})

#BOOST LOG
add_definitions(-DBOOST_ALL_DYN_LINK)
find_library(BOOST_LOG boost_log PATHS ${BOOST_LOG_LIB} NO_DEFAULT_PATH)
find_library(BOOST_LOG_SETUP boost_log_setup PATHS ${BOOST_LOG_LIB} NO_DEFAULT_PATH)
find_library(BOOST_DATE_TIME boost_date_time PATHS ${BOOST_LOG_LIB} NO_DEFAULT_PATH)
find_library(BOOST_SYSTEM boost_system PATHS ${BOOST_LOG_LIB} NO_DEFAULT_PATH)
find_library(BOOST_THREAD boost_thread PATHS ${BOOST_LOG_LIB} NO_DEFAULT_PATH)
find_library(BOOST_PO boost_program_options PATHS ${BOOST_LOG_LIB} NO_DEFAULT_PATH)
#BOOST LOG

#SuiteSparse
find_library(SS_AMD amd PATHS ${SS_CHOLMOD_LIB} NO_DEFAULT_PATH)
find_library(SS_CAMD camd PATHS ${SS_CHOLMOD_LIB} NO_DEFAULT_PATH)
find_library(SS_CCOLAMD ccolamd PATHS ${SS_CHOLMOD_LIB} NO_DEFAULT_PATH)
find_library(SS_CHOLMOD cholmod PATHS ${SS_CHOLMOD_LIB} NO_DEFAULT_PATH)
find_library(SS_COLAMD colamd PATHS ${SS_CHOLMOD_LIB} NO_DEFAULT_PATH)
find_library(SS_CONFIG suitesparseconfig PATHS ${SS_CHOLMOD_LIB} NO_DEFAULT_PATH)
find_library(SS_UMFPACK umfpack PATHS ${SS_CHOLMOD_LIB} NO_DEFAULT_PATH)
find_library(SS_QR spqr PATHS ${SS_CHOLMOD_LIB} NO_DEFAULT_PATH)
#SuiteSparse

find_library(TBB tbb PATHS ${TBB_LIB} NO_DEFAULT_PATH)
find_library(METIS metis PATHS ${METIS_LIB} NO_DEFAULT_PATH)
find_library(LAPACK lapack PATHS ${LAPACK_LIB} NO_DEFAULT_PATH)
find_library(BLAS openblas PATHS ${BLAS_LIB} NO_DEFAULT_PATH)

# TODO Currently testing/interface is done in main.cpp it should be moved to src directory where everything defined for lp
# Idea is to use this as linked library rather than header does it make sense!! lets see
include_directories(${CMAKE_SOURCE_DIR}/include)

# find_library(CONE_SOLVER coneSolver PATHS ${PROJECT_BINARY_DIR}/src)

# FIXME Move all find_library to src Cmake

add_executable(linearprogramming main.cpp)

target_link_libraries(${PROJECT_NAME} coneSolver
		      ${BOOST_LOG} ${BOOST_LOG_SETUP} ${BOOST_DATE_TIME} ${BOOST_PO}
		      ${BOOST_SYSTEM} ${BOOST_THREAD}
		      ${SS_QR}
		      ${SS_CHOLMOD}
		      ${SS_CONFIG}
		      ${METIS} ${SS_AMD} ${SS_CAMD} ${SS_COLAMD} ${SS_CCOLAMD} ${TBB}
		      ${BLAS} ${LAPACK} m rt
		      pthread hwloc gfortran)

install(TARGETS linearprogramming RUNTIME DESTINATION bin)


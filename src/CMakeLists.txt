
include_directories(${EIGEN_INC})
include_directories(${BOOST_INC})
include_directories(${SS_INC})

include_directories(${CMAKE_SOURCE_DIR}/include)

#BOOST LOG
add_definitions(-DBOOST_ALL_DYN_LINK)
find_library(BOOST_LOG boost_log PATHS ${BOOST_LOG_LIB} NO_DEFAULT_PATH)
find_library(BOOST_LOG_SETUP boost_log_setup PATHS ${BOOST_LOG_LIB} NO_DEFAULT_PATH)
find_library(BOOST_DATE_TIME boost_date_time PATHS ${BOOST_LOG_LIB} NO_DEFAULT_PATH)
find_library(BOOST_SYSTEM boost_system PATHS ${BOOST_LOG_LIB} NO_DEFAULT_PATH)
find_library(BOOST_THREAD boost_thread PATHS ${BOOST_LOG_LIB} NO_DEFAULT_PATH)
find_library(BOOST_PO boost_program_options PATHS ${BOOST_LOG_LIB} NO_DEFAULT_PATH)
#BOOST LOG

add_library(coneSolver STATIC LPFormatParser.cpp Problem.cpp Timer.cpp
                               ${BOOST_LOG} ${BOOST_LOG_SETUP} ${BOOST_DATE_TIME} ${BOOST_PO}
		               ${BOOST_SYSTEM} ${BOOST_THREAD} )
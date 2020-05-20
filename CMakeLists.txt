cmake_minimum_required(VERSION 3.17)
project(Efficient_Frontier)
message(STATUS "start running cmake...")

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_FLAGS "-std=c++11 -stdlib=libc++")
find_package(Boost 1.57.0 COMPONENTS system filesystem unit_test_framework REQUIRED)
find_package(Boost 1.57.0 COMPONENTS system filesystem unit_test_framework REQUIRED)
# Path to qunatlib library. On mac install with homebrew "brew install quantlib"
file(GLOB LIBRARIES "/usr/local/Cellar/quantlib/1.19/lib/*.a")
message("LIBRARIES = ${LIBRARIES}")
if(Boost_FOUND)

    message(STATUS "Boost_INCLUDE_DIRS: ${Boost_INCLUDE_DIRS}")
    message(STATUS "Boost_LIBRARIES: ${Boost_LIBRARIES}")
    message(STATUS "Boost_VERSION: ${Boost_VERSION}")

    include_directories(${Boost_INCLUDE_DIRS})
    ADD_DEFINITIONS(-DBOOST_TEST_DYN_LINK)
endif()

add_executable(Efficient_Frontier main.cpp)
if(Boost_FOUND)

    target_link_libraries(Efficient_Frontier ${Boost_LIBRARIES})
    target_link_libraries(Efficient_Frontier ${LIBRARIES})

endif()

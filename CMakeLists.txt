
cmake_minimum_required(VERSION 3.10)
project(forward_kinematics VERSION 1.0)

# specify the C++ standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# check if you are in the top-level project or not
if (NOT CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
    message(STATUS 
    "This project has a top-level one called 
    [${CMAKE_PROJECT_NAME}]")
else()
    message(STATUS 
    "This project is a top-level one")
endif()


add_subdirectory(kinetics)
add_subdirectory(utils)

# add_executable(${PROJECT_NAME} ../src/test_main.cpp)
# add_executable(${PROJECT_NAME} ../src/test_main_2.cpp)

add_executable(forward_kinematics ../src/test_main_4.cpp)
set (EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/bin)

# add_executable(test_kinematics ../src/test_kinematics.cpp)
# set (EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/bin)

target_include_directories(forward_kinematics 
# PUBLIC ${PROJECT_BINARY_DIR}
PUBLIC ${PROJECT_SOURCE_DIR}/include
)

# target_include_directories(test_kinematics 
# PUBLIC ${PROJECT_BINARY_DIR}
# PUBLIC ${PROJECT_SOURCE_DIR}/include
# )


link_directories(${PROJECT_SOURCE_DIR}/lib)

# target_link_libraries(${PROJECT_NAME} openGJKlib test)

target_link_libraries(forward_kinematics kinetics)
target_link_libraries(forward_kinematics utils)

# target_link_libraries(test_kinematics kinetics)
# target_link_libraries(test_kinematics utils)




# gdb调试
SET(CMAKE_BUILD_TYPE "Debug")  
SET(CMAKE_CXX_FLAGS_DEBUG "$ENV{CXXFLAGS} -O0 -Wall -g2 -ggdb")  
SET(CMAKE_CXX_FLAGS_RELEASE "$ENV{CXXFLAGS} -O3 -Wall") 


                           
                         








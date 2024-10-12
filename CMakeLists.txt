# You should always specify a range with the newest and oldest tested versions
# of CMake.
cmake_minimum_required(VERSION 3.22)

# set version
set(PROJECT_VERSION 0.0.0)
set(PROJECT_NAME EKF_DSSE)

# This is your project statement. You should always list languages;
project(${PROJECT_NAME} VERSION ${PROJECT_VERSION})

# specify the C++ standard
set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

include(FetchContent)
FetchContent_Declare(
  SuiteSparse
  GIT_REPOSITORY https://github.com/DrTimothyAldenDavis/SuiteSparse.git
  GIT_TAG v7.8.2)

FetchContent_Declare(
  json
  GIT_REPOSITORY https://github.com/nlohmann/json.git
  GIT_TAG v3.11.3)

FetchContent_MakeAvailable(SuiteSparse json)

add_subdirectory(libs)
add_subdirectory(state_estimator)

add_library(ekf_algo src/ekf.cpp)

add_compile_definitions("CS_COMPLEX")

# Create ALIAS targets.
add_library(${CMAKE_PROJECT_NAME}::ekf ALIAS ekf_algo)

target_include_directories(ekf_algo PUBLIC include)

target_link_libraries(ekf_algo PUBLIC SuiteSparse::CXSparse)

add_library(ekf_interface src/interface.cpp)

# Create ALIAS targets.
add_library(${CMAKE_PROJECT_NAME}::interface ALIAS ekf_interface)

set(DATA_INPUT_DIR "${CMAKE_SOURCE_DIR}/inputs")
set(DATA_OUTPUT_DIR "${CMAKE_SOURCE_DIR}/outputs")

# Make sure the output directory exists
file(MAKE_DIRECTORY "${DATA_OUTPUT_DIR}")

# Pass the directories to your program as compile definitions
target_compile_definitions(
  ekf_interface PRIVATE "-DDATA_INPUT_DIR=\"${DATA_INPUT_DIR}\""
                        "-DDATA_OUTPUT_DIR=\"${DATA_OUTPUT_DIR}\"")

target_include_directories(ekf_interface PUBLIC include)

target_link_libraries(ekf_interface PUBLIC nlohmann_json::nlohmann_json)

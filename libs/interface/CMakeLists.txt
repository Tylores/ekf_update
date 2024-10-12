add_library(interface src/interface.cpp)

# Set the data input directory
set(DATA_INPUT_DIR "${CMAKE_SOURCE_DIR}/inputs")

# Set the data output directory
set(DATA_OUTPUT_DIR "${CMAKE_SOURCE_DIR}/outputs")

# Make sure the output directory exists
file(MAKE_DIRECTORY "${DATA_OUTPUT_DIR}")

# Pass the directories to your program as compile definitions
target_compile_definitions(
  interface PRIVATE "-DDATA_INPUT_DIR=\"${DATA_INPUT_DIR}\""
                     "-DDATA_OUTPUT_DIR=\"${DATA_OUTPUT_DIR}\"")

target_include_directories(interface PUBLIC include)

target_link_libraries(interface PUBLIC nlohmann_json::nlohmann_json)

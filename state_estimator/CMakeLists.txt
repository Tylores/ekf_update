add_executable(state_estimator src/main.cpp)

add_compile_definitions("DEBUG")

set(DATA_INPUT_DIR "${CMAKE_SOURCE_DIR}/inputs")
set(DATA_OUTPUT_DIR "${CMAKE_SOURCE_DIR}/outputs")

# Make sure the output directory exists
file(MAKE_DIRECTORY "${DATA_OUTPUT_DIR}")

# Pass the directories to your program as compile definitions
target_compile_definitions(
  state_estimator PRIVATE "-DDATA_INPUT_DIR=\"${DATA_INPUT_DIR}\""
                        "-DDATA_OUTPUT_DIR=\"${DATA_OUTPUT_DIR}\"")

target_link_libraries(state_estimator 
  ${CMAKE_PROJECT_NAME}::ekf 
  ${CMAKE_PROJECT_NAME}::interface
)

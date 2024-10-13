add_executable(state_estimator src/main.cpp)

add_compile_definitions("DEBUG")

target_link_libraries(state_estimator ${CMAKE_PROJECT_NAME}::interface)

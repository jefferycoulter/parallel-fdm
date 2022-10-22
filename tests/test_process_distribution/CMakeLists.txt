project(process_distribution)

# set the output directory
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/tests/bin)

# specify source and header files
set(SRCS ${CMAKE_CURRENT_SOURCE_DIR}/test_process_distribution.c
    ${CMAKE_SOURCE_DIR}/src/fdm.c
    ${CMAKE_SOURCE_DIR}/src/io.c
    ${CMAKE_SOURCE_DIR}/src/parallel.c
)

# specify the executable to be built
add_executable(${PROJECT_NAME} ${SRCS})

# specify which directories to include
target_include_directories(${PROJECT_NAME} PUBLIC ${MPI_C_INCLUDE_DIRS} ${CMAKE_SOURCE_DIR}/inc)

# specify which libraries to link with
target_link_libraries(${PROJECT_NAME} PUBLIC ${MPI_C_LIBRARIES})
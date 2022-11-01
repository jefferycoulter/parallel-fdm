cmake_minimum_required(VERSION 3.16)

project(main)

add_compile_options(-g -Wall -Wextra)

find_package(MPI REQUIRED)
if(MPI_C_FOUND)
    message(STATUS "Found MPI")
endif(MPI_C_FOUND)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR} )

add_executable(${PROJECT_NAME} src/main.c
    inc/fdm.h
    inc/io.h
    inc/parallel.h
    inc/reaction.h
    inc/shape.h
    
    src/fdm.c
    src/io.c
    src/parallel.c
    src/reaction.c
)

target_include_directories(${PROJECT_NAME} PUBLIC ${MPI_C_INCLUDE_DIRS} inc)

target_link_libraries(${PROJECT_NAME} PUBLIC ${MPI_C_LIBRARIES})

add_subdirectory(tests)



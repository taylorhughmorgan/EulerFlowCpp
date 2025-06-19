cmake_minimum_required(VERSION 3.10)
project(MyCProject C)

# Add include directories
include_directories(${PROJECT_SOURCE_DIR}/include)
include_directories(${PROJECT_SOURCE_DIR}/lib)

# Add executable and source files
set(SOURCES
    main.c
    lib/rkf45.c
)

add_executable(main ${SOURCES})
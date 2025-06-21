cmake_minimum_required(VERSION 3.10)
project(MyCProject C)
find_package(GSL REQUIRED)
# Add include directories
include_directories(${PROJECT_SOURCE_DIR}/include)
include_directories(${PROJECT_SOURCE_DIR}/lib)

# Add executable and source files
set(SOURCES
    heat_eqn.c
)

add_executable(heat_eqn ${SOURCES})

# Link the math library
target_link_libraries(heat_eqn m)
target_link_libraries(heat_eqn GSL::gsl GSL::gslcblas)
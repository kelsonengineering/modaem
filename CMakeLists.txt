cmake_minimum_required(VERSION 3.20)

project(MODAEM
    VERSION 1.9.0
    DESCRIPTION "Analytic Element Model for groundwater flow"
    LANGUAGES Fortran
)

# Default to Release if not specified
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

# BLAS/LAPACK — prefer Apple's Accelerate on macOS
if(APPLE)
    set(BLA_VENDOR Apple)
endif()
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)

add_subdirectory(src)

install(TARGETS modaem
    RUNTIME DESTINATION bin
)

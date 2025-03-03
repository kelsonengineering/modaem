cmake_minimum_required(VERSION 3.12)

project(MODAEM LANGUAGES Fortran)

if(APPLE)
    set(BLA_VENDOR "Apple")
endif()

find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)

# Set Fortran standard
set(CMAKE_Fortran_STANDARD 90)
set(CMAKE_Fortran_STANDARD_REQUIRED ON)

# Set default build type to Debug if not specified
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Debug)
endif()

# Set compiler flags for different build types
set(CMAKE_Fortran_FLAGS_DEBUG "-g -O0")
set(CMAKE_Fortran_FLAGS_RELEASE "-O3")

# Create executable target
add_executable(modaem)

target_sources(modaem
    PRIVATE
        modaem.F90
        a_grid.F90
        a_obs.F90
        a_stdio.F90
        a_tr0.F90
        f_bwl.F90
        f_dipole.F90
        f_pond.F90
        f_well.F90
        i_areasink.F90
        i_dipole.F90
        i_linesink.F90
        i_pond.F90
        i_well.F90
        m_aem.F90
        m_aqu.F90
        m_as0.F90
        m_cw0.F90
        m_ext.F90
        m_hb0.F90
        m_in0.F90
        m_inq.F90
        m_ls0.F90
        m_ls1.F90
        m_ls2.F90
        m_ls3.F90
        m_pd0.F90
        m_wl0.F90
        m_wl1.F90
        u_constants.F90
        u_grid.F90
        u_io.F90
        u_math.F90
        u_matrix.F90
        u_polygon.F90
        ../../dependencies/polypack/polypack.f
    )

    target_link_libraries(modaem PRIVATE BLAS::BLAS LAPACK::LAPACK)

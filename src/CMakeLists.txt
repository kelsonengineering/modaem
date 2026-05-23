add_executable(modaem
    modaem.F90
    # Utility modules
    u_constants.F90
    u_domain.F90
    u_grid.F90
    u_io.F90
    u_math.F90
    u_matrix.F90
    u_polygon.F90
    # Functional elements
    f_bwl.F90
    f_dipole.F90
    f_linesink.F90
    f_pond.F90
    f_reference.F90
    f_well.F90
    # Interface modules
    i_areasink.F90
    i_dipole.F90
    i_linesink.F90
    i_pond.F90
    i_well.F90
    # Element modules
    f_aem.F90
    p_aqu.F90
    p_as0.F90
    p_cw0.F90
    a_ext.F90
    p_hb0.F90
    m_inq.F90
    p_ls0.F90
    p_ls1.F90
    p_ls2.F90
    p_ls3.F90
    p_pd0.F90
    p_wl0.F90
    p_wl1.F90
    # Application modules
    a_grid.F90
    a_obs.F90
    a_stdio.F90
    a_tr0.F90
    # Vendored dependencies
    dependencies/polypack/polypack.f
)

target_link_libraries(modaem PRIVATE BLAS::BLAS LAPACK::LAPACK)

# Compiler flags per build type
target_compile_options(modaem PRIVATE
    $<$<CONFIG:Debug>:-g -O0 -Wall -fcheck=all>
    $<$<CONFIG:Release>:-O3 -ftree-vectorize>
)

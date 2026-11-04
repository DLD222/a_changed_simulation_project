foreach(D IN LISTS WarpX_DIMS)
    warpx_set_suffix_dims(SD ${D})
    target_sources(lib_${SD}
      PRIVATE
        MomentFluidSolver.cpp
        MultiFluidContainer.cpp
        QdsmcParticleContainer.cpp
        WarpXFluidContainer.cpp
    )
endforeach()

# The 5-moment fluid reuses the header-only MomentSolver headers, whose
# internal includes are written relative to Source/MomentsSolver.
target_include_directories(lib_${SD} PRIVATE ${CMAKE_SOURCE_DIR}/Source/MomentsSolver)

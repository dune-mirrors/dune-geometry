dune_add_test(SOURCES test-affinegeometry.cc)

dune_add_test(SOURCES test-axisalignedcubegeometry.cc)

# whether the compiler is affected by
# https://gitlab.dune-project.org/core/dune-geometry/issues/19,
# a.k.a. https://gcc.gnu.org/bugzilla/show_bug.cgi?id=87288
set(ISSUE19_UNAFFECTED_GXX ON)
set(issue19_affected_configs)
function(check_compiler_issue19 config)
  string(TOUPPER "${config}" CONFIG)
  separate_arguments(opts UNIX_COMMAND "${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_${CONFIG}}")
  execute_process(
    COMMAND ${CMAKE_CXX_COMPILER} ${opts} -Q --help=optimizers
    OUTPUT_VARIABLE opt_settings
    ERROR_VARIABLE ignored # silence cc1 complaint about c++ options...
  )
  if(opt_settings MATCHES "-ftree-loop-vectorize[ \t]+\\[enabled\\]")
    set(ISSUE19_UNAFFECTED_GXX OFF PARENT_SCOPE)
    list(APPEND issue19_affected_configs "${config}")
    set(issue19_affected_configs "${issue19_affected_configs}" PARENT_SCOPE)
  endif()
endfunction(check_compiler_issue19)
if((CMAKE_CXX_COMPILER_ID STREQUAL "GNU") AND
    (CMAKE_CXX_COMPILER_VERSION VERSION_GREATER_EQUAL 8.2) AND
    (CMAKE_CXX_COMPILER_VERSION VERSION_LESS 8.3))
  if(CMAKE_BUILD_TYPE)
    check_compiler_issue19("${CMAKE_BUILD_TYPE}")
  elseif(CMAKE_CONFIGURATION_TYPES)
    # multi-configuration generator, check each configuration and if any is
    # bad skip test for all of them.  dune_add_test does not really allow to
    # make the guard configuration-dependent.
    foreach(conf IN LISTS CMAKE_CONFIGURATION_TYPES)
      check_compiler_issue19("${conf}")
    endforeach(conf)
  else()
    check_compiler_issue19("")
  endif()
  if(NOT ISSUE19_UNAFFECTED_GXX)
    if(issue19_affected_configs)
      string(REPLACE ";" ", " issue19_affected_configs "${issue19_affected_configs}")
      set(issue19_affected_configs
        " (Affected configurations: ${issue19_affected_configs})")
    endif()
    message(WARNING "You are using gcc 8.2 with optimization options "
      "including -ftree-loop-vectorize (possibly implicated by -O3).  This is "
      "known to miscompile, see "
      "https://gitlab.dune-project.org/core/dune-geometry/issues/19 and "
      "https://gcc.gnu.org/bugzilla/show_bug.cgi?id=87288${issue19_affected_configs}")
    message(WARNING "Skipping test-cornerstoragerefwrap due to known-bad "
      "compiler/optimization combination")
  endif()
endif()
dune_add_test(SOURCES test-cornerstoragerefwrap.cc
              CMAKE_GUARD ISSUE19_UNAFFECTED_GXX)

dune_add_test(SOURCES test-fromvertexcount.cc)

dune_add_test(SOURCES test-referenceelements.cc)

dune_add_test(SOURCES test-quadrature.cc)
add_dune_blas_lapack_flags(test-quadrature)

dune_add_test(SOURCES test-multilineargeometry.cc)

dune_add_test(SOURCES test-nonetype.cc)

dune_add_test(SOURCES test-refinement.cc)

dune_add_test(SOURCES test-constexpr-geometrytype.cc)

dune_add_test(SOURCES test-geometrytype-id.cc)

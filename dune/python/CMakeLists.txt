# if Python bindings are enabled, include necessary sub directories.
if (DUNE_ENABLE_PYTHONBINDINGS)
  add_subdirectory(geometry)
endif()
add_subdirectory(test)

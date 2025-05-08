# Common library sources
set(COMMON_SOURCES
    main.cpp
    octree.cpp
    philox_rng.cpp
)

# Create static library
add_library(nbody_common STATIC ${COMMON_SOURCES})

# Include directories
target_include_directories(nbody_common PUBLIC 
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CXX_IMPLICIT_INCLUDE_DIRECTORIES}
)

# Link against C++ standard library
target_link_libraries(nbody_common PUBLIC stdc++) 
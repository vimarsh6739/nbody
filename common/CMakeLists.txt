# Common library sources
set(COMMON_SOURCES
    octree.cpp
    philox_rng.cpp
)

# Create static library
add_library(nbody_common STATIC ${COMMON_SOURCES})


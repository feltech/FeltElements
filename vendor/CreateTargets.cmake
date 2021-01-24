include(ExternalProject)

###################################
# OpenVolumeMesh

# Download and install OpenVolumeMesh.
set(OpenVolumeMesh_DIR ${CMAKE_CURRENT_LIST_DIR}/OpenVolumeMesh)
ExternalProject_Add(
	OpenVolumeMeshGit
	GIT_REPOSITORY "https://www.graphics.rwth-aachen.de:9000/OpenVolumeMesh/OpenVolumeMesh.git"
	PREFIX "${CMAKE_CURRENT_BINARY_DIR}/OpenVolumeMesh"
	CMAKE_ARGS
	-DCMAKE_INSTALL_PREFIX:PATH=${OpenVolumeMesh_DIR}
	-DBUILD_SHARED_LIBS:BOOL=OFF
	BUILD_BYPRODUCTS # Required for ninja
)
# Create OpenVolumeMesh interface target.
add_library(EXT::OpenVolumeMesh INTERFACE IMPORTED)
set_property(
	TARGET EXT::OpenVolumeMesh
	PROPERTY INTERFACE_COMPILE_DEFINITIONS
	INCLUDE_TEMPLATES=1  # Required to avoid compile errors
)
set_property(
	TARGET EXT::OpenVolumeMesh
	PROPERTY INTERFACE_INCLUDE_DIRECTORIES ${OpenVolumeMesh_DIR}/include)
set_property(
	TARGET EXT::OpenVolumeMesh PROPERTY INTERFACE_LINK_DIRECTORIES ${OpenVolumeMesh_DIR}/lib)
set_property(
	TARGET EXT::OpenVolumeMesh PROPERTY INTERFACE_LINK_LIBRARIES OpenVolumeMesh)

###################################
# Fastor

# Create Fastor interface target.
add_library(EXT::Fastor INTERFACE IMPORTED)
set_property(
	TARGET EXT::Fastor
	PROPERTY INTERFACE_INCLUDE_DIRECTORIES ${CMAKE_CURRENT_LIST_DIR}/Fastor)

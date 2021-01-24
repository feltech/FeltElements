# Create static library
add_library(feltelements STATIC)
# Link OpenVolumeMesh
add_dependencies(feltelements OpenVolumeMeshGit)
# Add libraries to link.
# Note: this mix of libs requires `compiler.libcxx=libstdc++11` in ~/.conan/profiles/default
target_link_libraries(
	feltelements
	PUBLIC
	EXT::OpenVolumeMesh
	EXT::Fastor
	CONAN_PKG::eigen
	CONAN_PKG::spdlog
	PRIVATE
	project_warnings
	CONAN_PKG::fmt
	CONAN_PKG::boost
)
# Library include directories
target_include_directories(feltelements PUBLIC .)
# Put lib under a "lib" directory.
set_target_properties(feltelements PROPERTIES LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
# Add source files to target.
file(GLOB_RECURSE _src_lib LIST_DIRECTORIES false *.cpp)
target_sources(feltelements PRIVATE ${_src_lib})

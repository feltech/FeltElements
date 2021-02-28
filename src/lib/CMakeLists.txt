# Create static library
add_library(feltelements.lib STATIC)
set_target_properties(feltelements.lib PROPERTIES OUTPUT_NAME "feltelements")
find_package(OpenMP REQUIRED)
# Add libraries to link.
# Note: this mix of libs requires `compiler.libcxx=libstdc++11` in ~/.conan/profiles/default
target_link_libraries(
	feltelements.lib
	PUBLIC
	EXT::OpenVolumeMesh
	EXT::Fastor
	CONAN_PKG::eigen
	OpenMP::OpenMP_CXX
	PRIVATE
	project_warnings
	Boost::headers
	CONAN_PKG::fmt
	CONAN_PKG::spdlog
)
# Library include directories
target_include_directories(feltelements.lib PUBLIC .)
# Put lib under a "lib" directory.
set_target_properties(
	feltelements.lib PROPERTIES ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
# Add source files to target.
file(GLOB_RECURSE _src_lib LIST_DIRECTORIES false *.cpp)
target_sources(feltelements.lib PRIVATE ${_src_lib})

# Set default spdlog logging level.
if (CMAKE_BUILD_TYPE STREQUAL "Release")
	target_compile_definitions(feltelements.lib PRIVATE SPDLOG_ACTIVE_LEVEL=SPDLOG_LEVEL_WARN)
else()
	# Set default spdlog logging level.
	target_compile_definitions(feltelements.lib PRIVATE SPDLOG_ACTIVE_LEVEL=SPDLOG_LEVEL_DEBUG)
endif()

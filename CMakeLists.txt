cmake_minimum_required(VERSION 3.15)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

project(FeltElements)
# Enable testing
set(BUILD_TESTING)
enable_testing()

###################################
# Setup reusable interface projects
add_library(project_warnings INTERFACE)
include(CompilerWarnings)
set_project_warnings(project_warnings)

###################################
# Fetch conan CMake integration

if(NOT EXISTS "${CMAKE_BINARY_DIR}/conan.cmake")
	message(STATUS "Downloading conan.cmake from https://github.com/conan-io/cmake-conan")
	file(DOWNLOAD "https://github.com/conan-io/cmake-conan/raw/master/conan.cmake"
		"${CMAKE_BINARY_DIR}/conan.cmake")
endif()
include(${CMAKE_BINARY_DIR}/conan.cmake)
# Ensure conan is installed
conan_check(REQUIRED)

###################################
# Initialise targets

# Add extra warnings.
add_compile_options(
	# Enable SSE etc
	-march=native)
# Global compiler definitions
add_compile_definitions(
	# Required for OpenVolumeMesh
	INCLUDE_TEMPLATES=1)

# Create static library
add_library(feltelements STATIC)
# Create test executable
add_executable(tests)

###################################
# Install C++ libraries.

conan_cmake_run(
	REQUIRES
	boost/1.75.0
	catch2/2.11.0
	eigen/3.3.4@conan/stable
	fmt/6.2.0
	spdlog/1.6.1

	OPTIONS
	boost:header_only=True

	BASIC_SETUP # Run conan_basic_setup
	PROFILE default PROFILE_AUTO build_type # Use default profile override with Debug/Release
	BUILD missing  # build packages with no binary available
	CMAKE_TARGETS  # Generate CMake targets to link to
)
conan_load_buildinfo()

include(ExternalProject)
# Download and install OpenVolumeMesh
set(OpenVolumeMesh_DIR ${CMAKE_CURRENT_SOURCE_DIR}/vendor/OpenVolumeMesh)
ExternalProject_Add(
	OpenVolumeMeshGit
	GIT_REPOSITORY "https://www.graphics.rwth-aachen.de:9000/OpenVolumeMesh/OpenVolumeMesh.git"
	PREFIX "${CMAKE_CURRENT_BINARY_DIR}/OpenVolumeMesh"
	CMAKE_ARGS -DCMAKE_INSTALL_PREFIX:PATH=${OpenVolumeMesh_DIR}
	BUILD_BYPRODUCTS # Required for ninja
)
# Create OpenVolumeMesh interface target
add_library(EXT::OpenVolumeMesh INTERFACE IMPORTED)
set_property(
	TARGET EXT::OpenVolumeMesh
	PROPERTY INTERFACE_INCLUDE_DIRECTORIES ${OpenVolumeMesh_DIR}/include)
set_property(
	TARGET EXT::OpenVolumeMesh PROPERTY INTERFACE_LINK_DIRECTORIES ${OpenVolumeMesh_DIR}/lib)
if ("${CMAKE_BUILD_TYPE}" STREQUAL "Debug")
	set_property(
		TARGET EXT::OpenVolumeMesh PROPERTY INTERFACE_LINK_LIBRARIES OpenVolumeMeshd)
else()
	set_property(
		TARGET EXT::OpenVolumeMesh PROPERTY INTERFACE_LINK_LIBRARIES OpenVolumeMesh)
endif()

# Create Fastor interface target
add_library(EXT::Fastor INTERFACE IMPORTED)
set_property(
	TARGET EXT::Fastor
	PROPERTY INTERFACE_INCLUDE_DIRECTORIES ${CMAKE_CURRENT_SOURCE_DIR}/vendor/Fastor)

###################################
# Configure lib

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
target_include_directories(feltelements PUBLIC src/lib)
# Put lib under a "lib" directory.
set_target_properties(feltelements PROPERTIES LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
# Add source files to target.
file(GLOB_RECURSE _src_lib LIST_DIRECTORIES false src/lib/*.cpp)
target_sources(feltelements PRIVATE ${_src_lib})

###################################
# Configure tests

include(${CONAN_LIB_DIRS_CATCH2}/cmake/Catch2/Catch.cmake)
# Put tests executable under a "test" directory, rather than "bin".
set_target_properties(tests PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/test")
# Add source files to tests target.
file(GLOB_RECURSE _src_test LIST_DIRECTORIES false src/test/*.cpp)
target_sources(tests PRIVATE ${_src_test})
# Add sanitizer and coverage build options.
target_compile_options(
	tests PRIVATE -fsanitize=address,leak,undefined -static-libasan --coverage)
target_link_options(
	tests PRIVATE -fsanitize=address,leak,undefined -static-libasan --coverage)
# Add libraries to link.
# Note: this mix of libs requires `compiler.libcxx=libstdc++11` in ~/.conan/profiles/default
target_link_libraries(
	tests
	PRIVATE
	project_warnings
	feltelements
	CONAN_PKG::fmt
	CONAN_PKG::boost
	CONAN_PKG::catch2)
get_target_property(_test_runtime_dir tests RUNTIME_OUTPUT_DIRECTORY)
# Add CTest tests
catch_discover_tests(tests WORKING_DIRECTORY "${_test_runtime_dir}")
# Copy test resources
set(_resources_dir ${CMAKE_CURRENT_SOURCE_DIR}/resources/test)
# Gather all resources required for tests (i.e. mesh files).
file(
	GLOB_RECURSE _resources LIST_DIRECTORIES false RELATIVE ${_resources_dir}
	CONFIGURE_DEPENDS ${_resources_dir}/*
)
message(STATUS "Found resources in ${_resources_dir}: ${_resources}")
file(MAKE_DIRECTORY ${_test_runtime_dir}/resources)
foreach (_resource ${_resources})
	configure_file(
		${_resources_dir}/${_resource} ${_test_runtime_dir}/resources/${_resource} COPYONLY)
endforeach()
# Make directory to store artefacts created in tests for further debugging.
file(MAKE_DIRECTORY ${_test_runtime_dir}/artefacts)


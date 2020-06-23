cmake_minimum_required(VERSION 3.15)
include(ExternalProject)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

project(FeltElements)
# Enable testing
set(BUILD_TESTING)
enable_testing()

###################################
# Fetch conan CMake integration

if(NOT EXISTS "${CMAKE_BINARY_DIR}/conan.cmake")
	message(STATUS "Downloading conan.cmake from https://github.com/conan-io/cmake-conan")
	file(DOWNLOAD "https://github.com/conan-io/cmake-conan/raw/v0.14/conan.cmake"
		"${CMAKE_BINARY_DIR}/conan.cmake")
endif()
include(${CMAKE_BINARY_DIR}/conan.cmake)
# Ensure conan is installed
conan_check(REQUIRED)

###################################
# Initialise targets

# Add extra warnings.
add_compile_options(
	# Extra warnings
	-Wall -Wextra -Wpedantic
	# Enable SSE etc
	-march=native)
# Glob lib sources
file(GLOB_RECURSE _src_lib src/lib/*.cpp)
# Global compiler definitions
add_compile_definitions(
	# Required for OpenVolumeMesh
	INCLUDE_TEMPLATES=1)

# Create test executable
add_executable(tests)

###################################
# Install C++ libraries from conan.

conan_cmake_run(
	CONANFILE conanfile.txt  # Install from conanfile
	BASIC_SETUP CMAKE_TARGETS # Run conan_basic_setup
	PROFILE default PROFILE_AUTO build_type # Use default profile override with Debug/Release
	BUILD missing  # build packages with no binary available
)
conan_load_buildinfo()

set(OpenVolumeMesh_DIR ${CMAKE_CURRENT_SOURCE_DIR}/vendor/OpenVolumeMesh)
ExternalProject_Add(
	OpenVolumeMeshGit
	GIT_REPOSITORY "https://www.graphics.rwth-aachen.de:9000/OpenVolumeMesh/OpenVolumeMesh.git"
	PREFIX "${CMAKE_CURRENT_BINARY_DIR}/OpenVolumeMesh"
	CMAKE_ARGS -DCMAKE_INSTALL_PREFIX:PATH=${OpenVolumeMesh_DIR}
	BUILD_BYPRODUCTS # Required for ninja
)

# Library include directories
include_directories(src/lib SYSTEM vendor/Fastor ${OpenVolumeMesh_DIR}/include ${CONAN_INCLUDE_DIRS})

###################################
# Configure tests

include(${CONAN_LIB_DIRS_CATCH2}/cmake/Catch2/Catch.cmake)
# Put tests executable under a "test" directory, rather than "bin".
set_target_properties(tests PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/test")
# Add source files to tests target.
file(GLOB _src_test LIST_DIRECTORIES false src/test/*.cpp)
target_sources(tests PRIVATE ${_src_test} ${_src_lib})
# Add libraries to link.
conan_target_link_libraries(tests)
# Link OpenVolumeMesh
add_dependencies(tests OpenVolumeMeshGit)
target_link_directories(tests PRIVATE ${OpenVolumeMesh_DIR}/lib)
target_link_libraries(tests OpenVolumeMeshd)
# Add sanitizer and coverage build options.
target_compile_options(
	tests PRIVATE -fsanitize=address -static-libasan --coverage)
target_link_options(
	tests PRIVATE -fsanitize=address -static-libasan --coverage)
get_target_property(_test_runtime_dir tests RUNTIME_OUTPUT_DIRECTORY)
# Add CTest tests
catch_discover_tests(tests WORKING_DIRECTORY "${_test_runtime_dir}")
# Copy test resources
set(_resources_dir ${CMAKE_CURRENT_SOURCE_DIR}/src/test/resources)
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


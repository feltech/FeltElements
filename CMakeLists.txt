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
# Initialise global target config

add_compile_options(
	# Enable SSE etc. Note: increases performance but reduces CPU compatibility.
	-march=native)

###################################
# Install/configure conan libraries

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
	NO_OUTPUT_DIRS
)
conan_load_buildinfo()

###################################
# Install/configure external libraries

include (vendor/CreateTargets.cmake)

###################################
# Configure main FeltElements lib

add_subdirectory(src/lib)

###################################
# Configure tests

add_subdirectory(src/test)

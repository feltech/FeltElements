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
	eigen/3.3.9
	fmt/7.1.3
	spdlog/1.8.2
	
	OPTIONS
	# Frustratingly, to avoid linking all static libs (apart from the one we want), we must
	# explicitly exclude them all one by one.
	boost:without_atomic=True
	boost:without_chrono=True
	boost:without_container=True
	boost:without_context=True
	boost:without_contract=True
	boost:without_coroutine=True
	boost:without_date_time=True
	boost:without_exception=True
	boost:without_fiber=True
	boost:without_filesystem=True
	boost:without_graph=True
	boost:without_graph_parallel=True
	boost:without_iostreams=True
	boost:without_json=True
	boost:without_locale=True
	boost:without_log=True
	boost:without_math=True
	boost:without_mpi=True
	boost:without_nowide=True
	boost:without_program_options=False
	boost:without_python=True
	boost:without_random=True
	boost:without_regex=True
	boost:without_serialization=True
	boost:without_stacktrace=True
	boost:without_system=True
	boost:without_test=True
	boost:without_thread=True
	boost:without_timer=True
	boost:without_type_erasure=True
	boost:without_wave=True

	BASIC_SETUP # Run conan_basic_setup
	PROFILE default PROFILE_AUTO build_type # Use default profile override with Debug/Release
	BUILD missing  # build packages with no binary available
	CMAKE_TARGETS  # Generate CMake targets to link to
	NO_OUTPUT_DIRS  # Don't screw with specifically set output dirs
)
conan_load_buildinfo()

###################################
# Install/configure external libraries

include (vendor/CreateTargets.cmake)

###################################
# Configure main FeltElements lib

add_subdirectory(src/lib)

###################################
# Configure app

add_subdirectory(src/app)

###################################
# Configure tests

add_subdirectory(src/test)

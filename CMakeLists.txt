cmake_minimum_required(VERSION 3.15)
include(${CMAKE_ROOT}/Modules/ExternalProject.cmake)
include(${CMAKE_ROOT}/Modules/FetchContent.cmake)

project(FeltElements)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_VERBOSE_MAKEFILE ON)
set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

# Add extra warnings.
add_compile_options(-Wall -Wextra -Wpedantic)
enable_testing()

###################################
# Glob application files
file(GLOB_RECURSE _src_lib src/lib/*.cpp)
add_executable(tests src/test/main.cpp ${_src_lib})
set_target_properties(tests PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/test)

###################################
# Install C++ libraries from conan.

# Fetch conan CMake integration
if(NOT EXISTS "${CMAKE_BINARY_DIR}/conan.cmake")
	message(STATUS "Downloading conan.cmake from https://github.com/conan-io/cmake-conan")
	file(DOWNLOAD "https://github.com/conan-io/cmake-conan/raw/v0.14/conan.cmake"
		"${CMAKE_BINARY_DIR}/conan.cmake")
endif()
include(${CMAKE_BINARY_DIR}/conan.cmake)

# Ensure conan is installed
conan_check(REQUIRED)

conan_cmake_run(
	CONANFILE conanfile.txt  # Install from conanfile
	BASIC_SETUP # Run conan_basic_setup
	BUILD missing
	PROFILE default PROFILE_AUTO build_type # Use default profile override with Debug/Release
	BUILD missing
)
include(${CMAKE_CURRENT_BINARY_DIR}/conanbuildinfo.cmake)
conan_basic_setup(TARGETS)

###################################
# Configure test target
add_compile_definitions(
	# Required for OpenVolumeMesh
	INCLUDE_TEMPLATES=1
)
include_directories(src/lib ${CONAN_INCLUDE_DIRS})
include(${CONAN_LIB_DIRS_CATCH2}/cmake/Catch2/Catch.cmake)
conan_target_link_libraries(tests)
catch_discover_tests(tests)
target_compile_options(
	tests PRIVATE -fsanitize=thread -fsanitize-address-use-after-scope --coverage)
target_link_options(
	tests PRIVATE -fsanitize=thread -fsanitize-address-use-after-scope --coverage)
add_test(tests tests)
set(_resources_dir ${CMAKE_CURRENT_SOURCE_DIR}/src/test/resources)
get_target_property(_test_runtime_dir tests RUNTIME_OUTPUT_DIRECTORY)
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

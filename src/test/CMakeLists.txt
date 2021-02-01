# Create test executable
add_executable(feltelements.tests)
set_target_properties(feltelements.tests PROPERTIES OUTPUT_NAME "tests")
include(${CONAN_LIB_DIRS_CATCH2}/cmake/Catch2/Catch.cmake)
# Put tests executable under a "test" directory, rather than "bin".
set_target_properties(
	feltelements.tests PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/test")
# Add source files to tests target.
file(GLOB_RECURSE _src_test LIST_DIRECTORIES false *.cpp)
target_sources(feltelements.tests PRIVATE ${_src_test})
# Add sanitizer and coverage build options.
target_compile_options(
	feltelements.tests PRIVATE -fsanitize=address,leak,undefined -static-libasan --coverage)
target_link_options(
	feltelements.tests PRIVATE -fsanitize=address,leak,undefined -static-libasan --coverage)
# Add libraries to link.
# Note: this mix of libs requires `compiler.libcxx=libstdc++11` in ~/.conan/profiles/default
target_link_libraries(
	feltelements.tests
	PRIVATE
	project_warnings
	feltelements.lib

	CONAN_PKG::boost
	CONAN_PKG::catch2
	CONAN_PKG::fmt
	CONAN_PKG::spdlog)

get_target_property(_test_runtime_dir feltelements.tests RUNTIME_OUTPUT_DIRECTORY)
# Add CTest tests
catch_discover_tests(feltelements.tests WORKING_DIRECTORY "${_test_runtime_dir}")
# Copy test resources
set(_resources_dir ${CMAKE_SOURCE_DIR}/resources/test)
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


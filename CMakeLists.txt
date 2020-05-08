cmake_minimum_required(VERSION 3.15)
include(${CMAKE_ROOT}/Modules/ExternalProject.cmake)
include(${CMAKE_ROOT}/Modules/FetchContent.cmake)

project(FeltElements)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_VERBOSE_MAKEFILE ON)
set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

# Add extra warnings.
add_compile_options(-Wall -Wextra -Wpedantic)
###################################
# Create application targets.
file(GLOB_RECURSE _src_lib src/lib/*.cpp)
add_executable(app src/app/main.cpp ${_src_lib})
add_executable(tests src/test/main.cpp ${_src_lib})
set_target_properties(tests PROPERTIES RUNTIME_OUTPUT_DIRECTORY test)

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
# Install Lua libararies from LuaDist.
message(STATUS "Configuring LuaDist packages")
# Find luadist binary.
find_package(LuaDist REQUIRED)
# Specify libraries to install via a "luadistfile.txt" (to mirror "conanfile.txt").
set(LUADIST_FILE "luadistfile.txt" CACHE PATH "File containing newline-separated list of LuaDist packages to install")
# Ensure CMake reconfigures when luadistfile.txt is altered.
set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS ${LUADIST_FILE})
# Read the file and clean it up.
file(READ ${LUADIST_FILE} _libs)
string(STRIP ${_libs} _libs)
# Transform newline delimited to semicolon delimited, i.e. a CMake list.
string(REGEX REPLACE "\n" ";" _libs "${_libs}")
# Set path where to put LuaDist cache (includes src files for debug builds)
set(_luadist_cache_root ${CMAKE_CURRENT_SOURCE_DIR}/.luadist/${CMAKE_BUILD_TYPE})
# Set path for C++ libaries
set(RUNTIME_LIBRARY_DIR "lib" CACHE PATH "Location to install libraries relative to the build tree")
# Set path for Lua libaries (lua and shared libs).
set(RUNTIME_LUA_LIBRARY_DIR "${RUNTIME_LIBRARY_DIR}/lua" CACHE PATH
	"Location to install Lua libraries relative to the build tree")
set(_runtime_library_dir_abs "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/${RUNTIME_LIBRARY_DIR}")
set(_runtime_lua_library_dir_abs "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/${RUNTIME_LUA_LIBRARY_DIR}")
message(STATUS "LuaDist installing ${_libs}")
# TODO: uninstall extraneous packages + reinstall if variables change.
# Do the install.
execute_process(
	COMMAND ${LUADIST_BIN} ${_luadist_cache_root} install luajit ${_libs}
	-verbose -binary=false -source=true
	# Keep src when done - useful for debug builds, so we can breakpoint in the library source.
	-debug
	# Directory to install Lua libraries (pure Lua and compiled).
	-DINSTALL_LMOD=${_runtime_lua_library_dir_abs}
	-DINSTALL_CMOD=${_runtime_lua_library_dir_abs}
	# Directory to install C++ libraries.
	-DINSTALL_LIB=${_runtime_library_dir_abs}
	# Ensure C++ libraries can be found by Lua compiled modules at link time (e.g. lsqlite3.so -> libsqlite.so).
	-DCMAKE_PREFIX_PATH=${_runtime_library_dir_abs}
	# For C++ modules that depend on other C++ modules (e.g. lsqlite.so -> libsqlite.so).
	-DCMAKE_INSTALL_RPATH=${RUNTIME_LIBRARY_DIR}
	# Mirror build type of this project (i.e. Debug / Release).
	-DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
	# Store any error code.
	RESULT_VARIABLE _luadist_failed
)

# TODO: annoyingly, LuaDist exits with 1 if there is nothing to install.
if (_luadist_failed AND NOT _luadist_failed EQUAL 1)
	message(FATAL_ERROR "LuaDist failed with: ${_luadist_failed}")
endif()

###################################
# Configure application target.
find_library(LIBLUA lua PATHS ${_runtime_library_dir_abs})
# Add ad-hoc 3rd party C++ dependencies.
include_directories(src/lib vendor ${_luadist_cache_root}/include ${CONAN_INCLUDE_DIRS})
# Link installed libraries to main target.
target_link_libraries(app ${LIBLUA})
conan_target_link_libraries(app)
# Add a compile definition to the location of lua libraries in the binary directory.
add_compile_definitions(
	LUA_LIBRARY_DIR=\"${RUNTIME_LUA_LIBRARY_DIR}\"
	INCLUDE_TEMPLATES=1
)
###################################
# Configure test target
enable_testing()
include(${CONAN_LIB_DIRS_CATCH2}/cmake/Catch2/Catch.cmake)
conan_target_link_libraries(tests)
catch_discover_tests(tests)
target_compile_options(
	tests PRIVATE -fsanitize=thread -fsanitize-address-use-after-scope)
target_link_options(
	tests PRIVATE -fsanitize=thread -fsanitize-address-use-after-scope)
set(_resources_dir ${CMAKE_CURRENT_SOURCE_DIR}/resources)
get_target_property(_test_runtime_dir tests RUNTIME_OUTPUT_DIRECTORY)
file(
	GLOB_RECURSE _resources LIST_DIRECTORIES false RELATIVE ${_resources_dir}
	CONFIGURE_DEPENDS resources/*
)
message(STATUS "Found resources: ${_resources}")
file(MAKE_DIRECTORY ${_test_runtime_dir}/resources)
foreach (_resource ${_resources})
	configure_file(
		${_resources_dir}/${_resource} ${_test_runtime_dir}/resources/${_resource} COPYONLY)
endforeach()

###################################
# Check for moonscript files.
set(SCRIPTS_ROOT "${CMAKE_CURRENT_SOURCE_DIR}/src" CACHE PATH
	"Root directory to search for Lua/Moonscript files")
file(
	GLOB_RECURSE _scripts LIST_DIRECTORIES false RELATIVE ${SCRIPTS_ROOT}
	CONFIGURE_DEPENDS ${SCRIPTS_ROOT}/*.moon
)
# Ensure we can transpile the moonscript.
if (_scripts)
	message(STATUS "Configuring Moonscript transpilation")
	find_program(MOON moonc REQUIRED)
	if ("${MOON}" STREQUAL "MOON-NOTFOUND")
		message(FATAL_ERROR "moonscript binary (moonc) not found and is required to transpile")
	endif()
	# Add transpilation step to moonscripts.
	foreach (_script ${_scripts})
		get_filename_component(_dir ${_script} DIRECTORY)
		get_filename_component(_basename ${_script} NAME_WE)
		set(_script_path "${SCRIPTS_ROOT}/${_script}")
		# Mirror the source directory structure in the transpiled destination.
		set(_transpiled_path "${_runtime_lua_library_dir_abs}/${_dir}/${_basename}.lua")
		file(MAKE_DIRECTORY ${_runtime_lua_library_dir_abs}/${_dir})
		list(APPEND _transpiled_paths ${_transpiled_path})
		add_custom_command(
			OUTPUT ${_transpiled_path} DEPENDS ${_script_path}
			COMMAND ${MOON} -o ${_transpiled_path} ${_script_path}
			COMMENT "Transpiling ${_script} to ${_transpiled_path}"
			WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
		)
	endforeach ()
	# Gather individual moonscript transpilation commands under a target.
	add_custom_target(moonscripts DEPENDS ${_transpiled_paths})
	# Ensure moonscripts are compiled, if required, before building/running the application.
	add_dependencies(app moonscripts)
endif()

###################################
# Handle normal Lua scripts
file(
	GLOB_RECURSE _scripts LIST_DIRECTORIES false RELATIVE ${SCRIPTS_ROOT}
	CONFIGURE_DEPENDS ${SCRIPTS_ROOT}/*.lua
)
#message(STATUS "Found Lua scripts: ${_scripts}")
foreach (_script ${_scripts})
	configure_file(${SCRIPTS_ROOT}/${_script} ${_runtime_lua_library_dir_abs}/${_script} COPYONLY)
endforeach()


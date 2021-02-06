add_executable(feltelements.app main.cpp)
set_target_properties(feltelements.app PROPERTIES OUTPUT_NAME "feltelements")
# Put executable under a "app" directory, rather than "bin".
set_target_properties(feltelements.app PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/app")
target_link_libraries(
	feltelements.app
	PRIVATE
	project_warnings
	feltelements.lib

	CONAN_PKG::boost
	CONAN_PKG::fmt
	CONAN_PKG::spdlog
)

if (CMAKE_BUILD_TYPE STREQUAL "Release")
	# Enable interprocedural (whole-program) optimisation.
	include(CheckIPOSupported)
	check_ipo_supported(RESULT _ipo_supported OUTPUT error)
	if(_ipo_supported )
		message(STATUS "IPO / LTO enabled")
		set_property(TARGET feltelements.app PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)
	else()
		message(STATUS "IPO / LTO not supported: <${error}>")
	endif()
endif()
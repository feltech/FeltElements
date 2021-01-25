add_executable(feltelements.app main.cpp)
set_target_properties(feltelements.app PROPERTIES OUTPUT_NAME "feltelements")
# Put executable under a "app" directory, rather than "bin".
set_target_properties(feltelements.app PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/app")
target_link_libraries(
	feltelements.app
	PRIVATE
	feltelements.lib
	# CONAN_PKG::boost_program_options  # See top-level CMakeLists.txt.
	CONAN_PKG::fmt
)

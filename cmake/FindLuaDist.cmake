find_program(
	LUADIST_BIN luadist
	PATHS ../luadist
	PATH_SUFFIXES _install/bin
	DOC "LuaDist binary for installing Lua libraries"
)
INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(LUADIST DEFAULT_MSG LUADIST_BIN)
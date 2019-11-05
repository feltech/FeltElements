#include <iostream>

#define SOL_EXCEPTIONS_SAFE_PROPAGATION 1
#define SOL_ALL_SAFETIES_ON 1
#include <sol/sol.hpp>

#define LUA_LOCAL_PATH LUA_LIBRARY_DIR "/?.lua;" LUA_LIBRARY_DIR "/?"
#define LUA_LOCAL_CPATH LUA_LIBRARY_DIR "/?.so;" LUA_LIBRARY_DIR "/?"

int main(int, char *[])
{
	setenv("LUA_PATH", LUA_LOCAL_PATH, 1);
	setenv("LUA_CPATH", LUA_LOCAL_CPATH, 1);

	sol::state lua{};
	// open some common libraries
	lua.open_libraries(
		sol::lib::base, sol::lib::package, sol::lib::string, sol::lib::io, sol::lib::os,
		sol::lib::math, sol::lib::coroutine, sol::lib::table, sol::lib::debug);

	lua.safe_script_file(LUA_LIBRARY_DIR "/app/main.lua");

	return 0;
}

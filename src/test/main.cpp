#include <spdlog/sinks/stdout_color_sinks.h>
#include <spdlog/spdlog.h>

#include "util/Format.hpp"
#define CATCH_CONFIG_MAIN
#define CATCH_CONFIG_CONSOLE_WIDTH 200
#include <catch2/catch.hpp>

struct MyListener : Catch::TestEventListenerBase {
	using TestEventListenerBase::TestEventListenerBase; // inherit constructor
	void testRunStarting(Catch::TestRunInfo const&) override {
		spdlog::set_level(spdlog::level::debug);
		spdlog::set_default_logger(spdlog::stderr_color_mt("stderr"));
	}
};
CATCH_REGISTER_LISTENER( MyListener )

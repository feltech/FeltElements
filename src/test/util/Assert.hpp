#pragma once
#include <FeltElements/Typedefs.hpp>

namespace FeltElements::Test
{
constexpr Scalar epsilon = 0.00005;

auto equal = [](auto const & a, auto const & b) {
	using Tensor::Func::all_of;
	using Tensor::Func::abs;
	return all_of(abs(a - b) < epsilon);
};

template <class Tensor>
void check_equal(Tensor const & in, std::string_view const desc, Tensor const & expected)
{
	INFO(desc)
	INFO(in)
	CHECK(equal(in, expected));
}

template <class In, class Expected>
void check_equal(
	In const & in,
	std::string_view const desc_in,
	Expected const & expected,
	std::string_view const desc_expected)
{
	INFO(desc_in)
	INFO(in)
	INFO("==")
	INFO(desc_expected)
	INFO(expected)
	if constexpr (std::is_base_of_v<Eigen::MatrixBase<In>, In>)
	{
		CHECK(in.isApprox(expected, epsilon));
	}
	else
	{
		CHECK(equal(in, expected));
	}
}
}  // namespace FeltElements::Test
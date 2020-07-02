#pragma once

#include <fmt/format.h>
#include <fmt/ostream.h>

#include <boost/algorithm/string/join.hpp>
#include <eigen3/Eigen/Core>
#include <range/v3/view/subrange.hpp>
#include <range/v3/view/transform.hpp>

template <class Matrix>
using ostream_if_eigen = std::enable_if_t<!Matrix::IsVectorAtCompileTime, std::ostream &>;
template <class Matrix>
ostream_if_eigen<Matrix> operator<<(std::ostream & os, Matrix value)
{
	using namespace FeltElements::Tensor;
	using boost::algorithm::join;
	// std::to_string has non-configurable precision of too many decimal places.
	auto const constexpr to_string = [](auto const f) { return fmt::format("{:f}", f); };

	std::vector<std::string> is{};
	is.reserve(static_cast<size_t>(value.rows()));
	for (Eigen::Index i = 0; i < value.rows(); i++)
	{
		using Vec = Eigen::Matrix<typename Matrix::Scalar, Matrix::ColsAtCompileTime, 1>;
		Vec const & vec = value.row(i);
		is.push_back(join(
			ranges::subrange{vec.data(), vec.data() + vec.rows()} |
				ranges::views::transform(to_string),
			", "));
	}
	os << join(is, ",\n\t") << "\n";
	return os;
}

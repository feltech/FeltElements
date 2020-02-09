#include <ostream>
#include <unsupported/Eigen/CXX11/Tensor>
#include <range/v3/view/transform.hpp>
#include <boost/algorithm/string/join.hpp>

namespace Eigen
{

using Tensor1d = Tensor<double, 1>;

auto begin(Tensor1d & m)
{
	return m.data();
}

auto end(Tensor1d & m)
{
	return m.data() + m.size();
}

auto begin(Tensor1d const & m)
{
	return m.data();
}

auto end(Tensor1d const & m)
{
	return m.data() + m.size();
}

} // namespace Eigen

std::ostream& operator << (
	std::ostream& os, FeltElements::Tetrahedron::ElasticityTensor const& value)
{
	using Index = FeltElements::Tetrahedron::ElasticityTensor::Index;
	using boost::algorithm::join;
	constexpr auto dim = [](auto i) {
	  return Eigen::internal::get<
	      2, FeltElements::Tetrahedron::ElasticityTensor::Dimensions::Base>::value;
	};
	auto to_string = [](auto const f) {
		std::stringstream ss;
		ss << f;
		return ss.str();
	};

	std::array<std::string, dim(0)> is{};
	for (Index i = 0; i < dim(0); i++)
	{
		std::array<std::string, dim(1)> js{};
		for (Index j = 0; j < dim(1); j++)
		{
			std::array<std::string, dim(2)> ks{};
			for (Index k = 0; k < dim(2); k++)
			{
				Eigen::Tensor1d const & vec = value.chip(i, 0).chip(j, 0).chip(k, 0);
				ks[k] += "{" + join(ranges::views::transform(vec, to_string), ", ") + "}";
			}
			js[j] += "{" + join(ks, ", ") + "}";
		}
		is[i] += "\t" + join(js, ",\n\t");
	}
	std::string str = "{\n" + join(is, "\n}, {\n") + "\n}";
	os << str << std::endl;
	return os;
}
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

namespace
{
template <class Tensor, std::size_t idx>
const Eigen::Index size = Eigen::internal::get<idx, typename Tensor::Dimensions::Base>::value;

template <class Tensor, std::size_t idx>
constexpr auto dimension()
{
	return Eigen::internal::get<idx, typename Tensor::Dimensions::Base>::value;
};

template <class Tensor, std::size_t... dim>
constexpr auto dimensions(std::index_sequence<dim...>)
{
	return std::array<Eigen::Index, Tensor::NumIndices>{size<Tensor, dim>...};
};

template <class Tensor>
constexpr auto dimensions()
{
	return dimensions<Tensor>(std::make_index_sequence<Tensor::NumIndices>{});
}

// std::to_string has non-configurable precision of too many decimal places.
auto const to_string = [](auto const f) {
	std::stringstream ss;
	ss << f;
	return ss.str();
};

template <class Tensor, Eigen::Index N>
using ostream_if = std::enable_if_t<Tensor::NumIndices == N, std::ostream &>;
}

auto equal = [](auto const & a, auto const & b)
{
  Eigen::Tensor<bool, 0> comparison = ((a - b).abs() < 0.00001).all();
  return comparison(0);
};

template <class Tensor>
ostream_if<Tensor, 4> operator<< (std::ostream& os, Tensor const& value)
{
	using Index = typename Tensor::Index;
	using boost::algorithm::join;

	constexpr auto dims{dimensions<Tensor>()};

	std::array<std::string, dims[0]> is{};
	std::array<std::string, dims[1]> js{};
	std::array<std::string, dims[2]> ks{};
	for (Index i = 0; i < dims[0]; i++)
	{
		for (Index j = 0; j < dims[1]; j++)
		{
			for (Index k = 0; k < dims[2]; k++)
			{
				Eigen::Tensor1d const & vec = value.chip(i, 0).chip(j, 0).chip(k, 0);
				ks[k] = join(vec | ranges::views::transform(to_string), ", ");
			}
			js[j] = join(ks, "}, {");
		}
		is[i] = join(js, "}},\n\t{{");
	}
	os << "{\n\t{{" << join(is, "}}\n}, {\n\t{{") << "}}\n}\n";
	return os;
}

template <class Tensor>
ostream_if<Tensor, 2> operator<< (std::ostream& os, Tensor const& value)
{
	using Index = typename Tensor::Index;
	using boost::algorithm::join;

	constexpr auto dims{dimensions<Tensor>()};

	std::array<std::string, dims[0]> is{};
	for (Index i = 0; i < dims[0]; i++)
	{
		Eigen::Tensor1d const & vec = value.chip(i, 0);
		is[i] = join(vec | ranges::views::transform(to_string), ", ");
	}
	os << "{\n\t{" << join(is, "},\n\t{") << "}\n}";
	return os;
}

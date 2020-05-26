#include <ostream>

#include <unsupported/Eigen/CXX11/Tensor>
#include <range/v3/view/transform.hpp>
#include <range/v3/action/transform.hpp>
#include <range/v3/view/join.hpp>
#include <range/v3/to_container.hpp>
#include <boost/algorithm/string/join.hpp>

#include <FeltElements/TetGenIO.hpp>
#include <FeltElements/Tetrahedron.hpp>
#include <FeltElements/Attributes.hpp>

namespace
{
template <class Tensor, std::size_t idx>
const Eigen::Index size = Eigen::internal::get<idx, typename Tensor::Dimensions::Base>::value;

template <class Tensor, std::size_t... dim>
constexpr auto dimensions(std::index_sequence<dim...>)
{
	return std::array<Eigen::Index, Tensor::NumIndices>{size<Tensor, dim>...};
}

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

auto const to_int = [](auto const f) {
	return int{f};
};

template <class Tensor, Eigen::Index N>
using ostream_if = std::enable_if_t<Tensor::NumIndices == N, std::ostream &>;

template <class Tensor>
using InitList = typename Eigen::internal::Initializer<Tensor, Tensor::NumIndices>::InitList;
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
				Eigen::Tensor<FeltElements::Scalar, 1> const & vec =
					value.chip(i, 0).chip(j, 0).chip(k, 0);
				ks[k] = join(
					ranges::subrange{vec.data(), vec.data() + vec.size()} |
						ranges::views::transform(to_string), ", ");
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
		Eigen::Tensor<FeltElements::Scalar, 1> const & vec = value.chip(i, 0);
		is[i] = join(
			ranges::subrange{vec.data(), vec.data() + vec.size()} |
				ranges::views::transform(to_string), ", ");
	}
	os << "{\n\t{" << join(is, "},\n\t{") << "}\n}";
	return os;
}

inline std::ostream& operator<< (
	std::ostream& os, std::vector<OpenVolumeMesh::VertexHandle> const& vtxhs)
{
	using ranges::views::transform;
	using boost::algorithm::join;
	os << join(vtxhs | transform(to_int) | transform(to_string), ", ") << "\n";
	return os;
}

#include <OpenVolumeMesh/FileManager/FileManager.hh>
#include <catch2/catch.hpp>	 // Must come after `operator<<` definitions.

template <class Tensor>
void check_equal( Tensor const & in, std::string_view const desc, InitList<Tensor> const values)
{
	Tensor expected{};
	expected.setValues(values);
	INFO(desc)
	INFO(in)
	CHECK(equal(in, expected));
}

template <class Tensor>
void check_equal(
	Tensor const & in, std::string_view const desc_in,
	Tensor const & expected, std::string_view const desc_expected)
{
	INFO(desc_in)
	INFO(in)
	INFO("==")
	INFO(desc_expected)
	INFO(expected)
	CHECK(equal(in, expected));
}

inline auto load_tet(char const * const file_name)
{
	using namespace FeltElements;
	FeltElements::Mesh mesh;
	OpenVolumeMesh::IO::FileManager{}.readFile(file_name, mesh);
	auto const & vtxhs = Tetrahedron::vtxhs(mesh, 0);
	auto const & X = Tetrahedron::X(mesh, vtxhs);
	auto x = Tetrahedron::x(vtxhs, Tetrahedron::x(mesh));

	return std::tuple(X, x);
}

inline auto load_ovm_mesh(char const * const file_name)
{
	using namespace FeltElements;
	FeltElements::Mesh mesh;
	OpenVolumeMesh::IO::FileManager file_manager{};
	file_manager.readFile(file_name, mesh);
	return mesh;
}

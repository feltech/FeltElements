#include <FeltElements/Attributes.hpp>
#include <FeltElements/Derivatives.hpp>
#include <FeltElements/TetGenIO.hpp>
#include <boost/algorithm/string/join.hpp>
#include <ostream>
#include <range/v3/all.hpp>
#include <unsupported/Eigen/CXX11/Tensor>
#include <Fastor/Fastor.h>
#include <fmt/format.h>
#include <fmt/ostream.h>

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
	return fmt::format("{:f}", f);
};

auto const to_int = [](auto const f) {
	return int{f};
};

template <class Tensor, FeltElements::Tensor::Index N>
using ostream_if = std::enable_if_t<Tensor::dimension_t::value == N, std::ostream &>;

template <class Matrix>
using void_if_eigen = std::enable_if_t<!Matrix::IsVectorAtCompileTime, void>;
template <class Matrix>
using ostream_if_eigen = std::enable_if_t<!Matrix::IsVectorAtCompileTime, std::ostream &>;

constexpr FeltElements::Scalar epsilon = 0.00001;
}

auto equal = [](auto const & a, auto const & b)
{
  return Fastor::all_of(a - b < epsilon);
};

template <std::size_t dim0, std::size_t dim1, std::size_t dim2, std::size_t dim3>
std::ostream & operator<< (std::ostream& os, FeltElements::Tensor::Multi<dim0,dim1,dim2,dim3> value)
{
	using namespace FeltElements::Tensor;
	using boost::algorithm::join;

	std::array<std::string, value.dimension(0)> is{};
	std::array<std::string, value.dimension(1)> js{};
	std::array<std::string, value.dimension(2)> ks{};
	for (Index i = 0; i < value.dimension(0); i++)
	{
		for (Index j = 0; j < value.dimension(1); j++)
		{
			for (Index k = 0; k < value.dimension(2); k++)
			{
				using FeltElements::Tensor::Func::all;
				using Vec = Vector<value.dimension(3)>;
				Vec const & vec = value(i, j, k, all);

				ks[k] = join(
					ranges::subrange{vec.data(), vec.data() + Vec::size()} |
						ranges::views::transform(to_string), ", ");
			}
			js[j] = join(ks, "}, {");
		}
		is[i] = join(js, "}},\n\t{{");
	}
	os << "{\n\t{{" << join(is, "}}\n}, {\n\t{{") << "}}\n}\n";
	return os;
}

template <std::size_t dim0, std::size_t dim1>
std::ostream & operator<< (std::ostream& os, FeltElements::Tensor::Matrix<dim0,dim1> value)
{
	using namespace FeltElements::Tensor;
	using boost::algorithm::join;

	std::array<std::string, value.dimension(0)> is{};
	for (Index i = 0; i < value.dimension(0); i++)
	{
		using FeltElements::Tensor::Func::all;
		using Vec = Vector<value.dimension(1)>;
		Vec const & vec = value(i, all);
		is[i] = join(
			ranges::subrange{vec.data(), vec.data() + Vec::size()} |
				ranges::views::transform(to_string), ", ");
	}
	os << "{\n\t{" << join(is, "},\n\t{") << "}\n}";
	return os;
}

template <class Matrix>
ostream_if_eigen<Matrix> operator<< (std::ostream & os, Matrix value)
{
	using namespace FeltElements::Tensor;
	using boost::algorithm::join;

	std::vector<std::string> is{};
	is.reserve(value.rows());
	for (Eigen::Index i = 0; i < value.rows(); i++)
	{
		using Vec = Eigen::Matrix<typename Matrix::Scalar, Matrix::ColsAtCompileTime, 1>;
		Vec const & vec = value.row(i);
		is.push_back(join(
			ranges::subrange{vec.data(), vec.data() + vec.rows()} |
			ranges::views::transform(to_string), ", "));
	}
	os << join(is, ",\n\t") << "\n";
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

inline std::ostream& operator<< (
	std::ostream& os, FeltElements::Mesh::PointT const& vtx)
{
	using ranges::views::transform;
	using boost::algorithm::join;
	os << fmt::format("{{{}, {}, {}}}\n", vtx[0], vtx[1], vtx[2]);
	return os;
}

#include <OpenVolumeMesh/FileManager/FileManager.hh>
#include <catch2/catch.hpp>	 // Must come after `operator<<` definitions.

template <class Tensor>
void check_equal( Tensor const & in, std::string_view const desc, Tensor const & expected)
{
	INFO(desc)
	INFO(in)
	CHECK(equal(in, expected));
}


template <class In, class Expected>
void check_equal(
	In const & in, std::string_view const desc_in,
	Expected const & expected, std::string_view const desc_expected)
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

inline auto load_ovm_mesh(char const * const file_name)
{
	using namespace FeltElements;
	FeltElements::Mesh mesh;
	OpenVolumeMesh::IO::FileManager file_manager{};
	file_manager.readFile(file_name, mesh);
	return mesh;
}

inline auto load_tet(char const * const file_name)
{
	using namespace FeltElements;
	FeltElements::Mesh mesh = load_ovm_mesh(file_name);
	Element::Attribute::VertexHandles const attrib_vtxhs{mesh};
	Element::Attribute::MaterialPosition const attrib_X{mesh, attrib_vtxhs};
	Node::Attribute::SpatialPosition const attrib_x{mesh};
	Node::Positions X = attrib_X[0];  // De-constify.
	Node::Positions x = attrib_x.for_element(attrib_vtxhs[0]);

	return std::tuple(X, x);
}


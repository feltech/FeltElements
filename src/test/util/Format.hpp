#pragma once
#ifndef EIGEN_DEFAULT_IO_FORMAT
#define EIGEN_DEFAULT_IO_FORMAT Eigen::IOFormat(3, 0, ", ", ",\n", "", "", "", "")
#endif
#include <fmt/format.h>
#include <fmt/ostream.h>

#include <FeltElements/Typedefs.hpp>
#include <boost/algorithm/string/join.hpp>
#include <range/v3/view/subrange.hpp>
#include <range/v3/view/transform.hpp>

// std::to_string has non-configurable precision of too many decimal places.
auto const to_string = [](auto const f) { return fmt::format("{:f}", f); };

template <std::size_t dim0, std::size_t dim1, std::size_t dim2, std::size_t dim3>
std::ostream & operator<<(
	std::ostream & os, FeltElements::Tensor::Multi<dim0, dim1, dim2, dim3> value)
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
						ranges::views::transform(to_string),
					", ");
			}
			js[j] = join(ks, "}, {");
		}
		is[i] = join(js, "}},\n\t{{");
	}
	os << "{\n\t{{" << join(is, "}}\n}, {\n\t{{") << "}}\n}\n";
	return os;
}

template <std::size_t dim0, std::size_t dim1, std::size_t dim2>
std::ostream & operator<<(std::ostream & os, FeltElements::Tensor::Multi<dim0, dim1, dim2> value)
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
			using FeltElements::Tensor::Func::all;
			using Vec = Vector<value.dimension(2)>;
			Vec const & vec = value(i, j, all);

			js[j] = join(
				ranges::subrange{vec.data(), vec.data() + Vec::size()} |
					ranges::views::transform(to_string),
				", ");
		}
		is[i] = join(js, "},\n\t{");
	}
	os << "{\n\t{" << join(is, "}\n}, {\n\t{") << "}\n}\n";
	return os;
}

template <std::size_t dim0, std::size_t dim1>
std::ostream & operator<<(std::ostream & os, FeltElements::Tensor::Matrix<dim0, dim1> value)
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
				ranges::views::transform(to_string),
			", ");
	}
	os << "{\n\t{" << join(is, "},\n\t{") << "}\n}";
	return os;
}

inline std::ostream & operator<<(
	std::ostream & os, std::vector<OpenVolumeMesh::VertexHandle> const & vtxhs)
{
	using boost::algorithm::join;
	using ranges::views::transform;
	os << join(vtxhs | transform([](auto vtxh) { return vtxh.idx(); }) | transform(to_string), ", ")
	   << "\n";
	return os;
}

inline std::ostream & operator<<(std::ostream & os, FeltElements::Mesh::PointT const & vtx)
{
	using boost::algorithm::join;
	using ranges::views::transform;
	os << fmt::format("{{{}, {}, {}}}\n", vtx[0], vtx[1], vtx[2]);
	return os;
}

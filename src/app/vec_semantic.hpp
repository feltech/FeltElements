#pragma once
#include <fmt/format.h>

#include <boost/range/irange.hpp>
#include <boost/program_options.hpp>

namespace boost::program_options
{
template <FeltElements::Tensor::Index N>
class vec_semantic : public value_semantic
{
public:
	using value_type = FeltElements::Tensor::Vector<N>;
	using scalar_type = typename value_type::scalar_type;

	explicit vec_semantic(value_type * store_to) : m_store_to{store_to} {}

	[[nodiscard]] vec_semantic * value_name(const std::string_view & value_name)
	{
		m_value_name = value_name;
		return this;
	}

	[[nodiscard]] std::string name() const override
	{
		if (!m_value_name.empty())
			return m_value_name;

		return fmt::format("REAL^{}", N);
	}
	[[nodiscard]] unsigned int min_tokens() const override
	{
		return value_type::size();
	}
	[[nodiscard]] unsigned int max_tokens() const override
	{
		return value_type::size();
	}
	[[nodiscard]] bool is_composing() const override
	{
		return false;
	}
	vec_semantic * required()
	{
		m_is_required = true;
		return this;
	}
	[[nodiscard]] bool is_required() const override
	{
		return m_is_required;
	}
	void parse(boost::any & value_store, std::vector<std::string> const & new_tokens, bool)
		const override
	{
		value_type value;
		for (auto idx : boost::irange(value_type::size()))
			value(idx) = boost::lexical_cast<scalar_type>(new_tokens[idx]);
		value_store = value;
	}
	bool apply_default(boost::any &) const override
	{
		return false;
	}
	void notify(boost::any const & value_store) const override
	{
		*m_store_to = boost::any_cast<value_type>(value_store);
	}

private:
	value_type * m_store_to;
	bool m_is_required{false};
	std::string m_value_name;
};

template <FeltElements::Tensor::Index N>
vec_semantic<N> * vec_value(typename vec_semantic<N>::value_type * v)
{
	return new vec_semantic(v);
}
}  // namespace boost::program_options
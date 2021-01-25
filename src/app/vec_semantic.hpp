#include <boost/range/irange.hpp>

namespace boost::program_options
{
template <FeltElements::Tensor::Index N>
class vec_semantic : public value_semantic
{
public:
	using value_type = FeltElements::Tensor::Vector<N>;
	using scalar_type = typename value_type::scalar_type;

	explicit vec_semantic(value_type * store_to) : m_store_to{store_to} {}

	[[nodiscard]] std::string name() const override
	{
		std::string n = "REAL";
		for (auto idx : boost::irange(value_type::size()))
			n += " REAL";
		return n;
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
	void parse(boost::any & value_store, std::vector<std::string> const & new_tokens, bool utf8)
		const override
	{
		if (new_tokens.size() != value_type::size())
			throw std::logic_error{"ignored"};

		value_type value;
		for (auto idx : boost::irange(value_type::size()))
			value(idx) = boost::lexical_cast<scalar_type>(new_tokens[idx]);
		value_store = value;
	}
	bool apply_default(boost::any & value_store) const override
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
};

template <FeltElements::Tensor::Index N>
vec_semantic<N> * vec_value(typename vec_semantic<N>::value_type * v)
{
	return new vec_semantic(v);
}
}  // namespace boost::program_options
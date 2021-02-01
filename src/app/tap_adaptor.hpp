#pragma once

#include <boost/iterator/transform_iterator.hpp>
#include <boost/range/adaptor/argument_fwd.hpp>
#include <boost/range/iterator_range.hpp>

namespace boost::adaptors
{
// Provide a range for your return type
template <typename Value, typename Fn>
class tap_value
{
public:
	typedef const Value & result_type;

	explicit tap_value(Fn fn) : m_fn{fn} {}

	Value const & operator()(Value const & x) const
	{
		m_fn(x);
		return x;
	}

private:
	Fn m_fn;
};

template <typename Range, typename Fn>
class tap_range : public boost::iterator_range<boost::transform_iterator<
					  tap_value<typename boost::range_value<Range>::type, Fn>,
					  typename boost::range_iterator<Range>::type>>
{
private:
	typedef typename boost::range_value<Range>::type value_type;
	typedef typename boost::range_iterator<Range>::type iterator_base;
	typedef tap_value<value_type, Fn> tapper;
	typedef boost::transform_iterator<tapper, iterator_base> tap_iterator;
	typedef boost::iterator_range<tap_iterator> base_t;

public:
	tap_range(Range & rng, Fn fn)
		: base_t{
			  tap_iterator{boost::begin(rng), tapper{fn}},
			  tap_iterator{boost::end(rng), tapper{fn}}}
	{
	}
};

// Implement a holder class to hold the arguments required to construct the RangeAdaptor.
// The holder combines multiple parameters into one that can be passed as the right operand of
// operator|().

template <typename Fn>
class tap_holder : public boost::range_detail::holder<Fn>
{
public:
	explicit tap_holder(Fn&& fn) : boost::range_detail::holder<Fn>(fn) {}
	explicit tap_holder(Fn const & fn) : boost::range_detail::holder<Fn>(fn) {}
	void operator=(const tap_holder &) = delete;
};

// Define an instance of the holder with the name of the adaptor

static boost::range_detail::forwarder<tap_holder> tapped =
	boost::range_detail::forwarder<tap_holder>();

// Define operator|

template <typename SinglePassRange, typename Fn>
inline tap_range<SinglePassRange, Fn> operator|(SinglePassRange & rng, const tap_holder<Fn> & f)
{
	return tap_range<SinglePassRange, Fn>(rng, f.val);
}

template <typename SinglePassRange, typename Fn>
inline tap_range<const SinglePassRange, Fn> operator|(
	const SinglePassRange & rng, const tap_holder<Fn> & f)
{
	return tap_range<const SinglePassRange, Fn>(rng, f.val);
}
}  // namespace boost::adaptors
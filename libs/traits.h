/**
 * tlibs2
 * type traits library
 * @author Tobias Weber <tobias.weber@tum.de>, <tweber@ill.fr>
 * @date Nov-2014 -- 2018
 * @license GPLv3, see 'LICENSE' file
 * @desc Forked on 7-Nov-2018 from the privately and TUM-PhD-developed "tlibs" project (https://github.com/t-weber/tlibs).
 */

#ifndef __TLIBS2_TRAITS_H__
#define __TLIBS2_TRAITS_H__

#include <type_traits>
#include <vector>
#include <array>
#include <list>
#include <initializer_list>
#include <utility>


namespace tl2 {


// -----------------------------------------------------------------------------
template<class T>
struct remove_constref
{
	typedef typename std::remove_const<
		typename std::remove_reference<T>::type
			>::type type;
};

template<class T>
using remove_constref_t = typename remove_constref<T>::type;
// -----------------------------------------------------------------------------



// -----------------------------------------------------------------------------
/**
 * function call implementation
 */
template<class t_func,
	class t_arg = double, template<class ...> class t_cont = std::vector,
	std::size_t... idx>
t_arg _call_impl(t_func func, const t_cont<t_arg>& args,
	const std::integer_sequence<std::size_t, idx...>&)
{
	return func(args[idx]...);
}

/**
 * function call implementation (specialisation for std::array)
 */
template<class t_func, class t_arg, std::size_t... idx>
t_arg _call_impl(t_func func, const std::array<t_arg, sizeof...(idx)>& args,
	const std::integer_sequence<std::size_t, idx...>&)
{
	return func(args[idx]...);
}


/**
 * call a function with the args from an STL container
 */
template<std::size_t iNumArgs, class t_func,
	class t_arg = double, template<class ...> class t_cont = std::vector>
t_arg call(t_func func, const t_cont<t_arg>& args)
{
	using t_seq = std::make_integer_sequence<std::size_t, iNumArgs>;
	return _call_impl<t_func, t_arg, t_cont>(func, args, t_seq());
}

/**
 * call a function with the args from a std::array
 */
template<std::size_t iNumArgs, class t_func, class t_arg = double>
t_arg call(t_func func, const std::array<t_arg, iNumArgs>& args)
{
	using t_seq = std::make_integer_sequence<std::size_t, iNumArgs>;
	return _call_impl<t_func, t_arg>(func, args, t_seq());
}


// -----------------------------------------------------------------------------


template<typename t_arg, std::size_t ...idx>
using _t_fkt_vararg_impl = t_arg(*)(
	typename std::remove_reference<
		decltype(std::declval<t_arg*>()[idx])
	>::type...);

template<typename t_arg, std::size_t ...idx>
static _t_fkt_vararg_impl<t_arg, idx...>
_tstfkt_vararg(const std::integer_sequence<std::size_t, idx...>&)
{ return nullptr; /* not interested in return value, only its type */ }


/**
 * constructs a function type with 'iNumArgs' arguments: t_arg (*) (t_arg, t_arg, ...)
 */
template<typename t_arg, std::size_t iNumArgs>
using t_fkt_vararg = decltype(
	_tstfkt_vararg<t_arg>(
		std::make_integer_sequence<std::size_t, iNumArgs>()));
// -----------------------------------------------------------------------------

}
#endif

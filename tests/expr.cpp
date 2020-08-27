/**
 * expression test
 * @author Tobias Weber <tweber@ill.fr>
 * @date 28-mar-20
 * @license GPLv3, see 'LICENSE' file
 *
 * g++-10 -std=c++20 -I.. -o expr expr.cpp ../libs/log.cpp
 */

#define BOOST_TEST_MODULE Expr Test
#include <boost/test/included/unit_test.hpp>
namespace test = boost::unit_test;
namespace testtools = boost::test_tools;


#include "libs/expr.h"
#include "libs/str.h"
#include "libs/math20.h"


using t_types_real = std::tuple<double, float>;
using t_types_int = std::tuple<int, long>;


BOOST_AUTO_TEST_CASE_TEMPLATE(test_expr_real, t_real, t_types_real)
{
	static constexpr t_real eps = 1e-6;
	tl2::ExprParser<t_real> parser;

	t_real result = parser.parse("1 + 2*3");
	BOOST_TEST(tl2::equals<t_real>(result, 7, eps));

	result = parser.parse("4 + 5*6");
	BOOST_TEST(tl2::equals<t_real>(result, 34, eps));

	result = parser.parse(" - (sqrt(4)-5)^3 -  5/2 ");
	BOOST_TEST(tl2::equals<t_real>(result, 24.5, eps));

	result = parser.parse("-cos(sin(1.23*pi))^(-1.2 + 3.2)");
	BOOST_TEST(tl2::equals<t_real>(result, -0.6228, 1e-3));
}


BOOST_AUTO_TEST_CASE_TEMPLATE(test_expr_func, t_real, t_types_real)
{
	static constexpr t_real eps = 1e-6;

	auto tupres = tl2::eval_expr<std::string, t_real>("\t2 + \t2*3*4\n");

	BOOST_TEST(
		std::get<0>(tupres),
		tl2::equals<t_real>(std::get<1>(tupres), 14., eps));
}


BOOST_AUTO_TEST_CASE_TEMPLATE(test_expr_int, t_int, t_types_int)
{
	tl2::ExprParser<t_int> parser;

	t_int result = parser.parse("1 + 2*3");
	BOOST_TEST(result == 7);

	result = parser.parse("4 + 5*6");
	BOOST_TEST(result == 34);
}

/**
 * crystal matrix text
 * @author Tobias Weber <tweber@ill.fr>
 * @date feb-19
 * @license GPLv3, see 'LICENSE' file
 *
 * g++ -std=c++17 -fconcepts -o cryst cryst.cpp
 * g++ -std=c++17 -fconcepts -DUSE_LAPACK -I/usr/include/lapacke -I/usr/local/opt/lapack/include -L/usr/local/opt/lapack/lib -o cryst cryst.cpp -llapacke
 */

#define BOOST_TEST_MODULE Xtal Test
#include <boost/test/included/unit_test.hpp>
namespace test = boost::unit_test;
namespace testtools = boost::test_tools;


#include <iostream>
#include <vector>

#include "../../libs/_cxx20/math_algos.h"
using namespace m_ops;


using t_types = std::tuple<long double, double, float>;
BOOST_AUTO_TEST_CASE_TEMPLATE(test_xtal, t_real, t_types)
{
	using t_cplx = std::complex<t_real>;
	using t_vec = std::vector<t_real>;
	using t_mat = m::mat<t_real, std::vector>;
	using t_vec_cplx = std::vector<t_cplx>;
	using t_mat_cplx = m::mat<t_cplx, std::vector>;

	auto A = m::A_matrix<t_mat, t_real>(3., 4., 5., 80./180.*M_PI, 100./180.*M_PI, 60./180.*m::pi<t_real>);
	auto B = m::B_matrix<t_mat, t_real>(3., 4., 5., 80./180.*M_PI, 100./180.*M_PI, 60./180.*m::pi<t_real>);
	auto [B2, ok] = m::inv<t_mat>(A);
	B2 = 2.*m::pi<t_real> * m::trans<t_mat>(B2);

	std::cout << "A  = " << A << std::endl;
	std::cout << "B  = " << B << std::endl;
	std::cout << "B2 = " << B2 << std::endl;

	BOOST_TEST(ok);
	BOOST_TEST(m::equals(B, B2, std::numeric_limits<t_real>::epsilon()*10.));
	std::cout << std::endl;
}

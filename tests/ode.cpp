/**
 * test ode calculation
 * @author Tobias Weber <tweber@ill.fr>
 * @date 28-oct-20
 * @license GPLv3, see 'LICENSE' file
 *
 * g++ -std=c++20 -DUSE_LAPACK -I.. -I/usr/include/lapacke -Iext/lapacke/include -Lext/lapacke/lib -o ode ode.cpp -llapacke
 */

//#define BOOST_TEST_MODULE Eigenvector Test
//#include <boost/test/included/unit_test.hpp>
//namespace test = boost::unit_test;
//namespace testtools = boost::test_tools;

#include <iostream>
#include <vector>

#include "libs/math20.h"
using namespace tl2_ops;


//using t_types = std::tuple<double, float>;
//BOOST_AUTO_TEST_CASE_TEMPLATE(test_eig, t_real, t_types)

int main()
{
	using t_real = double;
	using t_cplx = std::complex<t_real>;
	using t_vec_cplx = tl2::vec<t_cplx, std::vector>;
	using t_mat_cplx = tl2::mat<t_cplx, std::vector>;


	auto coeff = tl2::create<t_mat_cplx>({
		0., 1.,
		1., 1. });

	auto f0 = tl2::create<t_vec_cplx>({1., 1.});
	t_cplx x0 = 0.;
	t_cplx n0 = 0.;

	std::cout << "coeff = " << coeff << std::endl;
	std::cout << "x0 = " << x0 << ", f0 = " << f0 << std::endl;
	std::cout << std::endl;


	for(t_cplx x=0; x.real()<10; x+=1)
	{
		std::cout << "x = " << x << std::endl;
		auto [ok, f] = tl2_la::odesys_const<t_mat_cplx, t_vec_cplx, t_cplx>(coeff, x, x0, f0);
		std::cout << "ok = " << std::boolalpha << ok << std::endl;

		for(std::size_t i=0; i<f.size(); ++i)
			std::cout << "f_" << i << " = " << f[i] << std::endl;
		std::cout << std::endl;
	}

	std::cout << std::endl;

	for(t_cplx n=0; n.real()<10; n+=1)
	{
		std::cout << "n = " << n << std::endl;
		auto [ok, f] = tl2_la::diffsys_const<t_mat_cplx, t_vec_cplx, t_cplx>(coeff, n, n0, f0);
		std::cout << "ok = " << std::boolalpha << ok << std::endl;

		for(std::size_t i=0; i<f.size(); ++i)
			std::cout << "f_" << i << " = " << f[i] << std::endl;
		std::cout << std::endl;
	}


	//BOOST_TEST(ok);
	//BOOST_TEST(params[0] == 3.9, testtools::tolerance(1e-3));
	//BOOST_TEST(params[1] == 1.036, testtools::tolerance(1e-3));

	return 0;
}

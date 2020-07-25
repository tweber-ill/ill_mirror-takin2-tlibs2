/**
 * test eigensystem calculation
 * @author Tobias Weber <tweber@ill.fr>
 * @date 25-jul-20
 * @license GPLv3, see 'LICENSE' file
 *
 * g++ -std=c++20 -DUSE_LAPACK -I/usr/include/lapacke -Iext/lapacke/include -Lext/lapacke/lib -o eig eig.cpp -llapacke
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
	using t_vec = tl2::vec<t_real, std::vector>;
	using t_mat = tl2::mat<t_real, std::vector>;
	using t_vec_cplx = tl2::vec<t_cplx, std::vector>;
	using t_mat_cplx = tl2::mat<t_cplx, std::vector>;

	auto mat = tl2::create<t_mat>({
		1.5, 0., 0., 
		0., 1.5, 0., 
		0., 0., 1.5 });

	auto [ok, evals_re, evals_im, evecs_re, evecs_im] =
		tl2_la::eigenvec<t_mat, t_vec, t_real>(mat, false, false, true);

	for(std::size_t i=0; i<evals_re.size(); ++i)
		std::cout << "Re(eval): " << evals_re[i] << std::endl;
	for(std::size_t i=0; i<evals_im.size(); ++i)
		std::cout << "Im(eval): " << evals_im[i] << std::endl;
	for(std::size_t i=0; i<evecs_re.size(); ++i)
		std::cout << "Re(evec): " << evecs_re[i] << std::endl;
	for(std::size_t i=0; i<evecs_im.size(); ++i)
		std::cout << "Im(evec): " << evecs_im[i] << std::endl;

	//BOOST_TEST(ok);
	//BOOST_TEST(params[0] == 3.9, testtools::tolerance(1e-3));
	//BOOST_TEST(params[1] == 1.036, testtools::tolerance(1e-3));
}

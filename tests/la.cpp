/**
 * lapack test
 * @author Tobias Weber <tweber@ill.fr>
 * @date feb-19
 * @license GPLv3, see 'LICENSE' file
 *
 * g++ -std=c++20 -DUSE_LAPACK -I/usr/include/lapacke -I/usr/local/opt/lapack/include -L/usr/local/opt/lapack/lib -o la la.cpp -llapacke
 */

#define BOOST_TEST_MODULE La1
#include <boost/test/included/unit_test.hpp>
namespace test = boost::unit_test;
namespace testtools = boost::test_tools;

#include <iostream>
#include <vector>

#include "libs/math20.h"
using namespace tl2_ops;

using t_types = std::tuple<double, float>;
BOOST_AUTO_TEST_CASE_TEMPLATE(test_equals, t_real, t_types)
{
	using t_cplx = std::complex<t_real>;
	using t_vec = std::vector<t_real>;
	using t_mat = tl2::mat<t_real, std::vector>;
	using t_vec_cplx = std::vector<t_cplx>;
	using t_mat_cplx = tl2::mat<t_cplx, std::vector>;


	auto M = tl2::create<t_mat>({1, 2, 3, 3, 2, 6, 4, 2, 4});
	auto Z = tl2::create<t_mat_cplx>({1, 2, 3, 3, 2, 6, 4, 2, 4});
	std::cout << "M = " << M << std::endl;
	std::cout << "Z = " << Z << std::endl;


	{
		auto [ok, Q, R] = tl2_la::qr<t_mat>(M);
		auto [ok2, P, L, U] = tl2_la::lu<t_mat>(M);

		auto QR = Q*R;
		auto PLU = P*L*U;

		std::cout << "\nok = " << std::boolalpha << ok << std::endl;
		std::cout << "Q = " << Q << std::endl;
		std::cout << "R = " << R << std::endl;
		std::cout << "QR = " << QR << std::endl;

		std::cout << "\nok2 = " << std::boolalpha << ok2 << std::endl;
		std::cout << "P = " << P << std::endl;
		std::cout << "L = " << L << std::endl;
		std::cout << "U = " << U << std::endl;
		std::cout << "PLU = " << PLU << std::endl;

		BOOST_TEST(ok);
		BOOST_TEST(ok2);
		BOOST_TEST(tl2::equals(QR, M, 1e-4));
		BOOST_TEST(tl2::equals(PLU, M, 1e-4));
	}


	{
		auto [ok, Q, R] = tl2_la::qr<t_mat_cplx>(Z);
		auto [ok2, P, L, U] = tl2_la::lu<t_mat_cplx>(Z);

		auto QR = Q*R;
		auto PLU = P*L*U;

		std::cout << "\nok = " << std::boolalpha << ok << std::endl;
		std::cout << "Q = " << Q << std::endl;
		std::cout << "R = " << R << std::endl;
		std::cout << "QR = " << QR << std::endl;

		std::cout << "\nok2 = " << std::boolalpha << ok2 << std::endl;
		std::cout << "P = " << P << std::endl;
		std::cout << "L = " << L << std::endl;
		std::cout << "U = " << U << std::endl;
		std::cout << "PLU = " << PLU << std::endl;

		BOOST_TEST(ok);
		BOOST_TEST(ok2);
		BOOST_TEST(tl2::equals(QR, Z, 1e-4));
		BOOST_TEST(tl2::equals(PLU, Z, 1e-4));
	}


	{
		auto [ok, evals, evecs] = tl2_la::eigenvec<t_mat_cplx, t_vec_cplx>(Z, 0, 0, 1);
		std::cout << "\nok = " << std::boolalpha << ok << std::endl;
		for(std::size_t i=0; i<evals.size(); ++i)
			std::cout << "eval: " << evals[i] << ", evec: " << evecs[i] << std::endl;


		auto [ok2, U, Vh, vals] = tl2_la::singval<t_mat_cplx>(Z);
		std::cout << "\nok = " << std::boolalpha << ok2 << std::endl;
		std::cout << "singvals: ";
		for(std::size_t i=0; i<vals.size(); ++i)
			std::cout << vals[i] << " ";
		std::cout << std::endl;
		std::cout << "U = " << U << "\nVh = " << Vh << std::endl;

		std::cout << "diag{vals} * UVh = " << U*tl2::diag<t_mat_cplx>(vals)*Vh << std::endl;


		auto [inva, ok3a] = tl2_la::pseudoinv<t_mat_cplx>(Z);
		auto [invb, ok3b] = tl2::inv<t_mat_cplx>(Z);
		std::cout << "\nok = " << std::boolalpha << ok3a << ", " << ok3b << std::endl;
		std::cout << "pseudoinv = " << inva << std::endl;
		std::cout << "      inv  = " << invb << std::endl;


		BOOST_TEST(ok);
		BOOST_TEST(ok2);
		BOOST_TEST(ok3a);
		BOOST_TEST(ok3b);
	}

	{
		auto [ok, evals_re, evals_im, evecs_re, evecs_im] = tl2_la::eigenvec<t_mat, t_vec>(M, 0, 0, 1);
		std::cout << "\nok = " << std::boolalpha << ok << std::endl;
		for(std::size_t i=0; i<evals_re.size(); ++i)
			std::cout << "eval: " << evals_re[i] << " + i*" << evals_im[i]
			<< ", evec: " << evecs_re[i] << " +i*" << evecs_im[i] << std::endl;


		auto [ok2, U, Vt, vals] = tl2_la::singval<t_mat>(M);
		std::cout << "\nok = " << std::boolalpha << ok2 << std::endl;
		std::cout << "singvals: ";
		for(std::size_t i=0; i<vals.size(); ++i)
			std::cout << vals[i] << " ";
		std::cout << std::endl;
		std::cout << "U = " << U << "\nVt = " << Vt << std::endl;

		std::cout << "diag{vals} * UVt = " << U*tl2::diag<t_mat>(vals)*Vt << std::endl;


		auto [inva, ok3a] = tl2_la::pseudoinv<t_mat>(M);
		auto [invb, ok3b] = tl2::inv<t_mat>(M);
		std::cout << "\nok = " << std::boolalpha << ok3a << ", " << ok3b << std::endl;
		std::cout << "pseudoinv = " << inva << std::endl;
		std::cout << "      inv  = " << invb << std::endl;


		BOOST_TEST(ok);
		BOOST_TEST(ok2);
		BOOST_TEST(ok3a);
		BOOST_TEST(ok3b);
	}
}

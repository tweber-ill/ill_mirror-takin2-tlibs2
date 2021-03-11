/**
 * math lib test
 * @author Tobias Weber <tweber@ill.fr>
 * @date mar-21
 * @license GPLv3, see 'LICENSE' file
 */

#define BOOST_TEST_MODULE La1
#include <boost/test/included/unit_test.hpp>
namespace test = boost::unit_test;
namespace testtools = boost::test_tools;

#include <iostream>
#include <vector>

#include "libs/maths.h"


using t_types = std::tuple<double, float>;
BOOST_AUTO_TEST_CASE_TEMPLATE(test_equals, t_real, t_types)
{
	using namespace tl2_ops;

	using t_cplx = std::complex<t_real>;
	using t_vec = tl2::vec<t_real, std::vector>;
	using t_mat = tl2::mat<t_real, std::vector>;
	using t_vec_cplx = tl2::vec<t_cplx, std::vector>;
	using t_mat_cplx = tl2::mat<t_cplx, std::vector>;


	{
		auto M = tl2::create<t_mat>({
			1., 2., 3.,
			3., 1., 4.,
			9., -4., 2
		});

		// test determinant
		t_real det = tl2::det(M);
		std::cout << "M = " << M << std::endl;
		std::cout << "|M| = " << det << std::endl;

		BOOST_TEST(tl2::equals<t_real>(det, 15., 1e-4));
	}


	{
		auto M = tl2::create<t_mat_cplx>({
			1., 2., 3.,
			3., 1., 4.,
			9., -4., 2
		});

		// test determinant
		t_cplx det = tl2::det(M);
		std::cout << "M = " << M << std::endl;
		std::cout << "|M| = " << det << std::endl;

		BOOST_TEST(tl2::equals<t_cplx>(det, 15., 1e-4));
	}
}

/**
 * math lib test
 * @author Tobias Weber <tweber@ill.fr>
 * @date mar-19
 * @license GPLv3, see 'LICENSE' file
 *
 * g++ -std=c++17 -fconcepts -o mat1 mat1.cpp
 * g++ -std=c++17 -fconcepts -DUSE_LAPACK -I/usr/include/lapacke -I/usr/local/opt/lapack/include -L/usr/local/opt/lapack/lib -o leastsq leastsq.cpp -llapacke
 */

#define BOOST_TEST_MODULE Mat1
#include <boost/test/included/unit_test.hpp>
namespace test = boost::unit_test;
namespace testtools = boost::test_tools;

#include <iostream>
#include <vector>

#include "../../libs/math_algos.h"
using namespace m_ops;


using t_types = std::tuple<double, float>;
BOOST_AUTO_TEST_CASE_TEMPLATE(test_equals, t_real, t_types)
{
	using t_cplx = std::complex<t_real>;
	using t_vec = std::vector<t_real>;
	using t_mat = m::mat<t_real, std::vector>;
	using t_vec_cplx = std::vector<t_cplx>;
	using t_mat_cplx = m::mat<t_cplx, std::vector>;


	std::cout << m::stoval<unsigned int>("123") << std::endl;

	std::vector vec1{{
		m::create<t_vec>({1, 2, 3, 4, 5, 6, 7, 8, 9, 10}),
		m::create<t_vec>({5, 5, 7, 9, 9.5, 10.5, 10.5, 12, 13.5, 14})
	}};
	std::vector vec2{{
		m::create<t_vec>({1, 2, 3, 4, 5, 6, 7, 8, 9, 10.5}),
		m::create<t_vec>({5, 5, 7, 9, 9.5, 10.5, 10.5, 12, 13.5, 14})
	}};
	std::vector vec3{{
		m::create<t_vec>({1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 10}),
		m::create<t_vec>({5, 5, 7, 9, 9.5, 10.5, 10.5, 12, 13.5, 14, 14})
	}};

	
	BOOST_TEST(m::equals_all(vec1, vec1, 1e-5));
	BOOST_TEST(m::equals_all(vec3, vec3, 1e-5));
	BOOST_TEST(!m::equals_all(vec1, vec2, 1e-5));
	BOOST_TEST(!m::equals_all(vec1, vec3, 1e-5));
}

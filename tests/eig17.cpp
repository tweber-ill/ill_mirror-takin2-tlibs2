/**
 * test eigensystem calculation
 * @author Tobias Weber <tweber@ill.fr>
 * @date 11-aug-20
 * @license GPLv3, see 'LICENSE' file
 *
 * g++ -std=c++17 -DUSE_LAPACK -I.. -I/usr/include/lapacke -Iext/lapacke/include -Lext/lapacke/lib -o eig17 eig17.cpp -llapacke
 */

#include <iostream>
#include <vector>

#include "libs/math17.h"


using t_real = double;
using t_cplx = std::complex<double>;
using t_mat = ublas::matrix<t_cplx>;
using t_vec = ublas::vector<t_cplx>;

int main()
{
	auto mat = tl2::make_mat<t_mat>({
		{1.5,  0.,   0.},
		{0.,   1.0,  0.},
		{0.,   0.,   0.5 }});

	std::cout << mat << std::endl;

	std::vector<t_vec> evecs;
	std::vector<t_real> evals;

	std::cout << std::boolalpha << tl2::eigenvecsel_herm(mat, evecs, evals, false) << std::endl;

	for(t_real val : evals)
		std::cout << "eval: " << val << std::endl;

	for(const t_vec& vec : evecs)
		std::cout << "evec: " << vec << std::endl;

	return 0;
}

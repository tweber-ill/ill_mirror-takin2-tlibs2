/**
 * runtime library using the tlibs2 math library
 * @author Tobias Weber <tweber@ill.fr>
 * @date 17-apr-20
 * @license see 'LICENSE' file
 * @desc Forked on 18/July/2020 from my privatly developed "matrix_calc" project (https://github.com/t-weber/matrix_calc).
 */

#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <cfloat>
#include <vector>

#include "math20.h"


using t_real =  double;
using t_int = int64_t;

using t_vec = tl2::vec<t_real, std::vector>;
using t_mat = tl2::mat<t_real, std::vector>;

static t_real g_eps = DBL_EPSILON;


extern "C" {

/**
 * set float epsilon
 */
void ext_set_eps(t_real eps)
{
	g_eps = eps;
}


/**
 * get float epsilon
 */
t_real ext_get_eps()
{
	return g_eps;
}


/**
 * tests equality of floating point numbers
 */
int ext_equals(t_real x, t_real y, t_real eps)
{
	return tl2::equals<t_real>(x, y, eps);
}


/**
 * removes a given row and column of a square matrix
 */
void ext_submat(const t_real* M, t_int N, t_real* M_new, t_int iremove, t_int jremove)
{
	t_vec mat(N*N, M);
	t_vec submat = tl2::flat_submat(mat, N, N, iremove, jremove);
	submat.to_array(M_new);
}


/**
 * calculates the determinant
 */
t_real ext_determinant(const t_real* M, t_int N)
{
	t_mat mat(N, N, M);
	return tl2::det<t_mat>(mat);
}


/**
 * inverted matrix
 */
t_int ext_inverse(const t_real* M, t_real* I, t_int N)
{
	t_real fullDet = ext_determinant(M, N);

	// fail if determinant is zero
	if(ext_equals(fullDet, 0., g_eps))
		return 0;

	t_mat mat(N, N, M);
	auto [inv, ok] = tl2::inv<t_mat>(mat);
	inv.to_array(I);

	return ok==true;
}


/**
 * matrix power
 */
t_int ext_power(const t_real* M, t_real* P, t_int N, t_int POW)
{
	t_mat mat(N, N, M);
	auto [result, ok] = tl2::pow(mat, POW);
	result.to_array(P);

	return ok==true;
}


/**
 * transposed matrix
 */
void ext_transpose(const t_real* M, t_real* T, t_int rows, t_int cols)
{
	for(t_int ctr=0; ctr<rows*cols; ++ctr)
	{
		t_int i = ctr/cols;
		t_int j = ctr%cols;
		T[j*rows + i] = M[i*cols + j];
	}
}

}

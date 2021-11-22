/**
 * tlibs2 -- (C-only) math library
 * @author Tobias Weber <tweber@ill.fr>
 * @date nov-2020
 * @note forked on 20-nov-2020 from the runtime library of my private "matrix_calc" project (https://github.com/t-weber/matrix_calc/blob/master/src/runtime.c).
 * @license GPLv2 or GPLv3, see 'LICENSE' file
 *
 * ----------------------------------------------------------------------------
 * tlibs
 * Copyright (C) 2017-2021  Tobias WEBER (Institut Laue-Langevin (ILL),
 *                          Grenoble, France).
 * Copyright (C) 2015-2017  Tobias WEBER (Technische Universitaet Muenchen
 *                          (TUM), Garching, Germany).
 * matrix_calc
 * Copyright (C) 2020       Tobias WEBER (privately developed).
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 or version 3 of the License.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 * ----------------------------------------------------------------------------
 */

#ifndef __TLIBS2_C_MATHLIB_H__
#define __TLIBS2_C_MATHLIB_H__


// ----------------------------------------------------------------------------
// linked list
// ----------------------------------------------------------------------------
struct tl2_list
{
	struct tl2_list *next;
	void *elem;
};

extern struct tl2_list* tl2_lst_create(void *elem);
extern struct tl2_list* tl2_lst_append(struct tl2_list *lst, void *elem);
extern void tl2_lst_remove(struct tl2_list *lst, void *elem);
extern void tl2_lst_free(struct tl2_list *lst);
// ----------------------------------------------------------------------------


// ----------------------------------------------------------------------------
// linalg functions
// ----------------------------------------------------------------------------
/**
 * set float epsilon
 */
extern void tl2_set_eps(double eps);

/**
 * get float epsilon
 */
extern double tl2_get_eps();

/**
 * tests equality of floating point numbers
 */
extern int tl2_flt_equals(double x, double y, double eps);

/**
 * set matrix elements to zero
 */
extern void tl2_mat_zero(double* M, int I, int J);

/**
 * set vector elements to zero
 */
extern void tl2_vec_zero(double* vec, int N);

/**
 * copy a vector
 */
extern void tl2_vec_cpy(double* dst, const double* src, int N);

/**
 * copy a matrix
 */
extern void tl2_mat_cpy(double* DST, const double* SRC, int I, int J);

/**
 * removes a given row and column of a square matrix
 */
extern void tl2_submat(const double* M, int N, double* M_new, int iremove, int jremove);

/**
 * calculates the determinant
 */
extern double tl2_determinant(const double* M, int N);

/**
 * inverted matrix
 */
extern int tl2_inverse(const double* M, double* I, int N);

/**
 * matrix-matrix product
 */
extern void tl2_matmat_mul(const double* M1, const double* M2, double *RES, int I, int J, int K);

/**
 * matrix-vector product
 */
extern void tl2_matvec_mul(const double* M, const double* v, double *res, int I, int J);

/**
 * transposed matrix
 */
extern void tl2_transpose(const double* M, double* T, int rows, int cols);

/**
 * vector inner product
 */
extern double tl2_inner(const double* v0, const double* v1, int N);

/**
 * vector outer product
 */
extern void tl2_outer(const double* v0, const double* v1, double *M, int N);

/**
 * 3-vector cross product
 */
extern void tl2_cross(const double* v0, const double* v1, double *res);

/**
 * vector length
 */
extern double tl2_vec_len(const double* vec, int N);

/**
 * vector addition
 */
extern void tl2_vec_add(const double* v0, const double* v1, double *res, int N);

/**
 * vector subtraction
 */
extern void tl2_vec_sub(const double* v0, const double* v1, double *res, int N);

/**
 * negative vector
 */
extern void tl2_vec_neg(const double* vec, double *res, int N);

/**
 * vector-scalar multiplication
 */
extern void tl2_vec_mul(const double* v, double s, double *res, int N);

/**
 * vector-scalar division
 */
extern void tl2_vec_div(const double* v, double s, double *res, int N);

/**
 * matrix addition
 */
extern void tl2_mat_add(const double* M0, const double* M1, double *RES, int I, int J);

/**
 * matrix subtraction
 */
extern void tl2_mat_sub(const double* M0, const double* M1, double *RES, int I, int J);

/**
 * negative matrix
 */
extern void tl2_mat_neg(const double* M, double *RES, int I, int J);

/**
 * matrix-scalar multiplication
 */
extern void tl2_mat_mul(const double* M, double s, double *RES, int I, int J);

/**
 * matrix-scalar division
 */
extern void tl2_mat_div(const double* M, double s, double *RES, int I, int J);

/**
 * mean vector
 */
extern void tl2_vec_mean(const struct tl2_list* veclist, const struct tl2_list* problist,
	double* mean, int N);

/**
 * covariance matrix
 */
extern void tl2_covariance(const struct tl2_list* veclist, const struct tl2_list* problist,
	double* COV, double* mean, int N);

/**
 * matrix trafo
 */
extern void tl2_mat_trafo(const double* M, const double* T, double* RES, int N, int ortho);

/**
 * resolution matrix
 */
extern void tl2_reso(const struct tl2_list* veclist, const struct tl2_list* problist,
	double* COV, double* RESO);
// ----------------------------------------------------------------------------


#endif

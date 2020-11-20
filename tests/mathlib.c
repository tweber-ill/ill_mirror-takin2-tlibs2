/**
 * tlibs2 -- (C-only) math library test
 * @author Tobias Weber <tweber@ill.fr>
 * @date nov-2020
 * @license GPLv2 or GPLv3, see 'LICENSE' file
 *
 * gcc -Wall -Wextra -std=c99 -o mathlib mathlib.c ../libs/mathlib.c
 */

#include "../libs/mathlib.h"
#include <stdio.h>
#include <stdlib.h>


int main()
{
	struct tl2_list* vecs = tl2_lst_create(0);
	struct tl2_list* probs = tl2_lst_create(0);
	vecs->elem = calloc(3, sizeof(double));
	probs->elem = calloc(1, sizeof(double));

	((double*)vecs->elem)[0] = 1.;
	((double*)vecs->elem)[1] = 2.;
	((double*)vecs->elem)[2] = 3.;
	*((double*)probs->elem) = 1.;


	struct tl2_list* vecs1 = tl2_lst_append(vecs, 0);
	struct tl2_list* probs1 = tl2_lst_append(probs, 0);
	vecs1->elem = calloc(3, sizeof(double));
	probs1->elem = calloc(1, sizeof(double));

	((double*)vecs1->elem)[0] = 3.;
	((double*)vecs1->elem)[1] = 2.;
	((double*)vecs1->elem)[2] = 1.;
	*((double*)probs1->elem) = 1.;


	struct tl2_list* vecs2 = tl2_lst_append(vecs, 0);
	struct tl2_list* probs2 = tl2_lst_append(probs, 0);
	vecs2->elem = calloc(3, sizeof(double));
	probs2->elem = calloc(1, sizeof(double));

	((double*)vecs2->elem)[0] = 5.;
	((double*)vecs2->elem)[1] = -2.;
	((double*)vecs2->elem)[2] = 0.5;
	*((double*)probs2->elem) = 0.5;


	double *mean = calloc(3, sizeof(double));
	tl2_vec_mean(vecs, probs, mean, 3);
	printf("mean: %f %f %f\n", mean[0], mean[1], mean[2]);
	free(mean);


	double *cov = calloc(3*3, sizeof(double));
	tl2_covariance(vecs, probs, cov, 3);
	printf("cov:\n\t%f %f %f\n\t%f %f %f\n\t%f %f %f\n",
		cov[0], cov[1], cov[2],
		cov[3], cov[4], cov[5],
		cov[6], cov[7], cov[8]);
	free(cov);


	tl2_lst_free(vecs);
	tl2_lst_free(probs);
	return 0;
}

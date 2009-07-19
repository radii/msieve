/*--------------------------------------------------------------------
This source distribution is placed in the public domain by its author,
Jason Papadopoulos. You may use it for any purpose, free of charge,
without having to notify anyone. I disclaim any responsibility for any
errors.

Optionally, please be nice and tell me if you find this source to be
useful. Again optionally, if you add to the functionality present here
please consider making those additions public too, so that others may 
benefit from your work.	
       				   --jasonp@boo.net 4/3/09
--------------------------------------------------------------------*/

#include <polyroot.h>

#define MAX_ITER 100

static uint32 find_one_root(dd_complex_t *poly, uint32 degree,
				dd_complex_t x, dd_complex_t *root,
				double eps) {
	uint32 i, j;
	dd_complex_t g, g2, gp, gm, h, sq, dx;
	dd_t abp, abm;
	double m = (double)degree;

	if (degree == 1) {
		root->r = dd_div_dd(dd_neg(poly[0].r), poly[1].r);
		root->i = dd_set_d(0.0);
		return 0;
	}

	for (i = 0; i < MAX_ITER; i++) {
		dd_complex_t b = poly[degree];
		dd_complex_t d = cplx_set_d(0.0, 0.0);
		dd_complex_t f = cplx_set_d(0.0, 0.0);
		dd_t err = cplx_norm(b);
		dd_t abx = cplx_norm(x);

		for (j = degree - 1; (int32)j >= 0; j--) {
			f = cplx_add(cplx_mul(x, f), d);
			d = cplx_add(cplx_mul(x, d), b);
			b = cplx_add(cplx_mul(x, b), poly[j]);
			err = dd_add_dd(cplx_norm(b),
					dd_mul_dd(abx, err));
		}

		err = dd_mul_d(err, eps);
		if (dd_cmp_dd(cplx_norm(b), err) <= 0) {
			*root = x;
			return 0;
		}

		g = cplx_div(d, b);
		g2 = cplx_mul(g, g);
		h = cplx_sub(g2, cplx_mul_d(cplx_div(f, b), 2.0));
		sq = cplx_sub(cplx_mul_d(h, m), g2);
		sq = cplx_sqrt(cplx_mul_d(sq, m - 1));

		gp = cplx_add(g, sq);
		abp = cplx_norm(gp);
		gm = cplx_sub(g, sq);
		abm = cplx_norm(gm);
		if (dd_cmp_dd(abp, abm) < 0) {
			gp = gm;
			abp = abm;
		}

		if (dd_cmp_d(abp, 0.0) == 0)
			return 1;

		dx = cplx_div(cplx_set_d(m, 0.0), gp);
		x = cplx_sub(x, dx);
	}

	return 2;
}

static void deflate_poly(dd_complex_t *poly, uint32 degree,
			dd_complex_t root) {
	uint32 i;
	dd_complex_t tmp_coeff = poly[degree];

	for (i = degree - 1; (int32)i >= 0; i--) {
		dd_complex_t t = poly[i];
		poly[i] = tmp_coeff;
		tmp_coeff = cplx_add(cplx_mul(root, tmp_coeff), t);
	}
}

uint32 find_poly_roots(dd_t *poly, uint32 degree, dd_complex_t *roots) {

	uint32 i, j, iter;
	dd_complex_t start_coeff[MAX_ROOTFINDER_DEGREE + 1];
	dd_complex_t coeff[MAX_ROOTFINDER_DEGREE + 1];
	double start_r, start_i;

	for (i = 0; i <= degree; i++) {
		start_coeff[i].r = poly[i];
		start_coeff[i].i = dd_set_d(0.0);
		coeff[i] = start_coeff[i];
	}
	i = degree;
	j = 0;
	start_r = 0.0513;
	start_i = 0.0918;

	while (i) {

		dd_complex_t next_root, polished_root;

		/* start the iteration from an initial point near
		   the origin, to increase the odds of finding the root
		   of coeff(x) that has the smallest magnitude. We
		   may have to try multiple starting guesses, since 
		   the derivatives of poly must not be zero at the
		   starting point */

		for (iter = 0; iter < MAX_ITER; iter++) {
			if (find_one_root(coeff, i, 
					cplx_set_d(start_r, start_i),
					&next_root, 1e-30) == 0) {
				break;
			}

			start_r += 0.0513;
			start_i += 0.0918;
		}
		if (iter == MAX_ITER)
			return 2;

		/* change roots with very small imaginary part to
		   be explicitly real roots */

		if (dd_cmp_dd(dd_fabs(next_root.i),
			      dd_mul_d(dd_fabs(next_root.r), 1e-30)) <= 0) {
			next_root.i = dd_set_d(0.0);
		}

		/* polish the root. With an initial tolerance close
		   to that of the machine epsilon, this step is only
		   necessary if there was some kind of numerical 
		   instability deflating the polynomial, and we expect
		   this to be very unlikely */

		deflate_poly(coeff, i--, next_root);
		if (find_one_root(start_coeff, degree,
				next_root, &polished_root, 1e-30) != 0)
			return 1;

		/* save the polished root, along with its complex
		   comjugate if the root is complex */

		roots[j++] = polished_root;
		if (dd_cmp_d(next_root.i, 0.0) != 0) {
			next_root.i = dd_neg(next_root.i);
			deflate_poly(coeff, i--, next_root);

			roots[j].r = polished_root.r;
			roots[j].i = dd_neg(polished_root.i);
			j++;
		}
	}

	return 0;
}

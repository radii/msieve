/*--------------------------------------------------------------------
This source distribution is placed in the public domain by its author,
Jason Papadopoulos. You may use it for any purpose, free of charge,
without having to notify anyone. I disclaim any responsibility for any
errors.

Optionally, please be nice and tell me if you find this source to be
useful. Again optionally, if you add to the functionality present here
please consider making those additions public too, so that others may 
benefit from your work.	

$Id$
--------------------------------------------------------------------*/

#include "common.h"

#define PHI 1.618033988749
#define COMP_PHI (2.0 - PHI)

/*-------------------------------------------------------------------------*/
static double evaluate(double *base, double dist, 
			double *search_dir, uint32 ndim,
			objective_func callback, void *extra)
{
	uint32 i;
	double curr_pt[MAX_VARS];

	for (i = 0; i < ndim; i++)
		curr_pt[i] = base[i] + dist * search_dir[i];

	return callback(curr_pt, extra);
}

/*-------------------------------------------------------------------------*/
static void bracket_min(double *base, double *search_dir,
			double *a_out, double *b_out, double *c_out,
			double *fb_out, uint32 ndim, 
			objective_func callback, void *extra)
{
	double fa, fb, fc, t;
	double a = *a_out;
	double b = *b_out;
	double c;

	fa = evaluate(base, a, search_dir, ndim, callback, extra);
	fb = evaluate(base, b, search_dir, ndim, callback, extra);
	if (fb > fa) {
		t = a; a = b; b = t;
		t = fa; fa = fb; fb = t;
	}
	c = b + PHI * (b - a);
	fc = evaluate(base, c, search_dir, ndim, callback, extra);

	while (fb > fc) {
		double r = (b - a) * (fb - fc);
		double q = (b - c) * (fb - fa);
		double u, max_u, fu;

		t = fabs(q - r);
		t = MAX(t, 1e-20);
		if (q - r < 0.0)
			t = -t;

		u = b - ((b - c) * q - (b - a) * r) / (2.0 * t);
		max_u = b + 100.0 * (c - b);

		if ((b - u) * (u - c) > 0.0) {
			fu = evaluate(base, u, search_dir, ndim, 
					callback, extra);

			if (fu < fc) {
				*a_out = b; *b_out = u; *c_out = c;
				*fb_out = fu;
				return;
			}
			else if (fu > fb) {
				*a_out = a; *b_out = b; *c_out = u;
				*fb_out = fb;
				return;
			}

			u = c + PHI * (c - b);
			fu = evaluate(base, u, search_dir, ndim,
					callback, extra);
		}
		else if ((c - u) * (u - max_u) > 0.0) {
			fu = evaluate(base, u, search_dir, ndim,
					callback, extra);

			if (fu < fc) {
				b = c; c = u; u = c + PHI * (c - b);
				fb = fc; 
				fc = fu; 
				fu = evaluate(base, u, search_dir, ndim,
						callback, extra);
			}
		}
		else if ((u - max_u) * (max_u - c) >= 0.0) {
			u = max_u;
			fu = evaluate(base, u, search_dir, ndim,
					callback, extra);
		}
		else {
			u = c + PHI * (c - b);
			fu = evaluate(base, u, search_dir, ndim,
					callback, extra);
		}

		a = b; b = c; c = u;
		fa = fb; fb = fc; fc = fu;
	}

	*a_out = a; *b_out = b; *c_out = c;
	*fb_out = fb;
}

/*-------------------------------------------------------------------------*/
static double minimize_line_core(double *base, double *search_dir,
			double a_in, double b_in, double c_in, 
			double fb_in, double tol, double *min_out, 
			uint32 ndim, int32 *status,
			objective_func callback, void *extra)
{
	int32 i;
	double a, b, x, w, v, u;
	double fx, fw, fv, fu;
	double e = 0.0;

	a = MIN(a_in, c_in);
	b = MAX(a_in, c_in);
	x = w = v = b_in;
	fx = fw = fv = fb_in;

	for (i = 0; i < 100; i++) {
		double xm = 0.5 * (a + b);
		double tol1 = tol * fabs(x) + 1e-20;
		double tol2 = 2.0 * tol1;
		double d = 0, etemp;

		if (fabs(x - xm) <= tol2 - 0.5 * (b - a)) {
			*min_out = x;
			return fx;
		}

		if (fabs(e) > tol1) {
			double r = (x - w) * (fx - fv);
			double q = (x - v) * (fx - fw);
			double p = (x - v) * q - (x - w) * r;
			
			q = 2.0 * (q - r);
			if (q > 0.0)
				p = -p;
			q = fabs(q);
			etemp = e;
			e = d;
			if (fabs(p) >= fabs(0.5 * q * etemp) ||
			    p <= q * (a - x) || p >= q * (b - x)) {
	
				e = b - x;
				if (x >= xm)
					e = a - x;				
				d = COMP_PHI * e;
			}
			else {
				d = p / q;
				u = x + d;
				if (u - a < tol2 || b - u < tol2)
					d = (xm >= x) ? tol1 : -tol1;
			}
		}
		else {
			e = (x >= xm) ? a - x : b - x;
			d = COMP_PHI * e;
		}

		if (fabs(d) < tol1)
			u = x + ((d >= 0) ? tol1 : -tol1);
		else
			u = x + d;
		fu = evaluate(base, u, search_dir, ndim, callback, extra);

		if (fu <= fx) {
			if (u >= x)
				a = x;
			else
				b = x;

			v = w; w = x; x = u;
			fv = fw; fw = fx; fx = fu;
		}
		else {
			if (u < x)
				a = u;
			else
				b = u;
			if (fu <= fw || w == x) {
				v = w; w = u;
				fv = fw; fw = fu;
			}
			else if (fu <= fv || v == x || v == w) {
				v = u;
				fv = fu;
			}
		}
	}

	printf("too many line iterations\n");
	*min_out = x;
	*status = 1;
	return fx;
}

/*-------------------------------------------------------------------------*/
static double minimize_line(double *base, double *search_dir, 
				uint32 ndim, int32 *status,
				objective_func callback,
				void *extra)
{
	uint32 i;
	double a = 0.0;
	double b = 1.0;
	double c, fb, min_dist, fmin_dist;

	bracket_min(base, search_dir, &a, &b, &c, 
			&fb, ndim, callback, extra);
	fmin_dist = minimize_line_core(base, search_dir, a, b, c, 
					fb, 1e-3, &min_dist, ndim,
					status, callback, extra);
	if (*status)
		return 0;

	for (i = 0; i < ndim; i++) {
		search_dir[i] *= min_dist;
		base[i] += search_dir[i];
	}
	return fmin_dist;
}

/*-------------------------------------------------------------------------*/
double minimize(double p[MAX_VARS], uint32 ndim, 
			double ftol, uint32 max_iter,
			objective_func callback, void *extra)
{
	uint32 i, j;
	double best_f, curr_f;
	double p_old[MAX_VARS], p_extrap[MAX_VARS], dir_extrap[MAX_VARS];
	double directions[MAX_VARS][MAX_VARS];
	int32 status = 0;

	memset(directions, 0, sizeof(directions));
	for (i = 0; i < ndim; i++) {
		p_old[i] = p[i];
		directions[i][i] = 1.0;
	}
	best_f = callback(p, extra);

	if (ndim == 1) {
		curr_f = minimize_line(p, directions[0], ndim,
					&status, callback, extra);
		if (status)
			return best_f;
		return curr_f;
	}

	for (i = 0; i < max_iter; i++) {

		int32 ibig = 0;
		double del = 0.0;
		double start_f = best_f;

		for (j = 0; j < ndim; j++) {
			curr_f = minimize_line(p, directions[j], ndim,
						&status, callback, extra);
			if (status)
				return best_f;

			if (best_f - curr_f > del) {
				ibig = j;
				del = best_f - curr_f;
			}
			best_f = curr_f;
		}

		if (2.0 * fabs(start_f - best_f) <= 
				ftol * (fabs(start_f) + fabs(best_f)) + 1e-20)
			return best_f;

		for (j = 0; j < ndim; j++) {
			p_extrap[j] = 2.0 * p[j] - p_old[j];
			dir_extrap[j] = p[j] - p_old[j];
			p_old[j] = p[j];
		}

		curr_f = callback(p_extrap, extra);
		if (curr_f < start_f) {

			double t1 = start_f - best_f - del;
			double t2 = start_f - curr_f;
			double t = 2.0 * (start_f - 2.0 * best_f + curr_f) *
					t1 * t1 - del * t2 * t2;

			if (t < 0.0) {
				best_f = minimize_line(p, dir_extrap, ndim,
							&status, callback, 
							extra);
				if (status)
					return best_f;

				for (j = 0; j < ndim; j++) {
					directions[ibig][j] = 
						directions[ndim-1][j];
					directions[ndim-1][j] = 
						dir_extrap[j];
				}
			}
		}
	}

	return best_f;
}

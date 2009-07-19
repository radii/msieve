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

#ifndef _DDCOMPLEX_H_
#define _DDCOMPLEX_H_

#include <util.h>
#include <dd.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
	dd_t r, i;
} dd_complex_t;

static INLINE dd_complex_t cplx_set(dd_t re, dd_t im) {
	dd_complex_t c;
	c.r = re;
	c.i = im;
	return c;
}

static INLINE dd_complex_t cplx_set_d(double re, double im) {
	return cplx_set(dd_set_d(re), dd_set_d(im));
}

static INLINE dd_complex_t cplx_add(dd_complex_t a, dd_complex_t b) {
	return cplx_set(dd_add_dd(a.r, b.r),
			dd_add_dd(a.i, b.i));
}

static INLINE dd_complex_t cplx_sub(dd_complex_t a, dd_complex_t b) {
	return cplx_set(dd_sub_dd(a.r, b.r),
			dd_sub_dd(a.i, b.i));
}

static INLINE dd_complex_t cplx_mul(dd_complex_t a, dd_complex_t b) {
	return cplx_set(dd_sub_dd(dd_mul_dd(a.r, b.r),
				  dd_mul_dd(a.i, b.i)),
			dd_add_dd(dd_mul_dd(a.r, b.i),
				  dd_mul_dd(a.i, b.r)));
}

static INLINE dd_complex_t cplx_mul_d(dd_complex_t a, double b) {
	return cplx_set(dd_mul_d(a.r, b),
			dd_mul_d(a.i, b));
}

static INLINE dd_t cplx_norm(dd_complex_t a) {
	dd_t xr = dd_fabs(a.r);
	dd_t xi = dd_fabs(a.i);
	if (dd_cmp_d(xr, 0.0) == 0) {
		return dd_set_dd(xi.hi, xi.lo);
	}
	else if (dd_cmp_d(xi, 0.0) == 0) {
		return dd_set_dd(xr.hi, xr.lo);
	}
	else if (dd_cmp_dd(xr, xi) > 0) {
		dd_t q = dd_div_dd(xi, xr);
		dd_t arg = dd_add_d(dd_mul_dd(q, q), 1.0);
		return dd_mul_dd(xr, dd_sqrt(arg));
	}
	else {
		dd_t q = dd_div_dd(xr, xi);
		dd_t arg = dd_add_d(dd_mul_dd(q, q), 1.0);
		return dd_mul_dd(xi, dd_sqrt(arg));
	}
}

static INLINE dd_complex_t cplx_div(dd_complex_t a, dd_complex_t b) {

	dd_complex_t ans;
	if (dd_cmp_dd(dd_fabs(b.r), dd_fabs(b.i)) >= 0) {
		dd_t q = dd_div_dd(b.i, b.r);
		dd_t den = dd_add_dd(b.r, dd_mul_dd(q, b.i));
		ans = cplx_set(dd_div_dd(dd_add_dd(a.r, 
						    dd_mul_dd(q, a.i)),
					  den),
				dd_div_dd(dd_sub_dd(a.i,
						    dd_mul_dd(q, a.r)),
					  den));
	}
	else {
		dd_t q = dd_div_dd(b.r, b.i);
		dd_t den = dd_add_dd(b.i, dd_mul_dd(q, b.r));
		ans = cplx_set(dd_div_dd(dd_add_dd(dd_mul_dd(q, a.r),
						    a.i), 
					  den),
				dd_div_dd(dd_sub_dd(dd_mul_dd(q, a.i),
						    a.r),
					  den));
	}

	return ans;
}

static INLINE dd_complex_t cplx_sqrt(dd_complex_t a) {

	dd_t xr, xi, w, r;
	dd_complex_t ans;

	if (dd_cmp_d(a.r, 0.0) == 0 &&
	    dd_cmp_d(a.i, 0.0) == 0)
		return a;
	
	xr = dd_fabs(a.r);
	xi = dd_fabs(a.i);
	if (dd_cmp_dd(xr, xi) >= 0) {
		r = dd_div_dd(xi, xr);
		w = dd_sqrt(dd_add_d(dd_mul_dd(r, r), 1.0));
		w = dd_mul_dpow2(dd_add_d(w, 1.0), 0.5);
		w = dd_mul_dd(dd_sqrt(xr), dd_sqrt(w));
	}
	else {
		r = dd_div_dd(xr, xi);
		w = dd_sqrt(dd_add_dd(dd_mul_dd(r, r), r));
		w = dd_mul_dpow2(dd_add_dd(w, r), 0.5);
		w = dd_mul_dd(dd_sqrt(xi), dd_sqrt(w));
	}

	if (dd_cmp_d(a.r, 0.0) >= 0) {
		ans.r = w;
		ans.i = dd_div_dd(a.i, dd_mul_dpow2(w, 2.0));
	}
	else {
		if (dd_cmp_d(a.i, 0.0) >= 0)
			ans.i = w;
		else
			ans.i = dd_set_dd(-w.hi, -w.lo);
		ans.r = dd_div_dd(a.i, dd_mul_d(ans.i, 2.0));
	}

	return ans;
}


#ifdef __cplusplus
}
#endif

#endif /* !_DDCOMPLEX_H_ */

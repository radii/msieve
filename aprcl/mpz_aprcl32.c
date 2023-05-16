/* Copyright 2011,2012,2013 David Cleaver
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

/*
 * v1.0 Posted to SourceForge on 2013/07/04
 *
 * v1.1 Posted to SourceForge on 2013/12/27
 * [The following fix was recommended by Dana Jacobsen and verified by Jon Grantham]
 *      - Bug fix: Removed unnecessary vl==0 check in mpz_extrastronglucas_prp
 * [The following improvements/fixes were recommended by Laurent Desnogues in 2013/08]
 *      - Speed improvement 1: Removed extraneous NormalizeJS calls in ARPCL
 *      - Speed improvement 2: Removed/consolidated calls to mpz_mod in APRCL
 *        (these improvements make the APRCL code about 1.5-2.2x faster)
 *      - Bug fix: Final test in APRCL routine is now correct
 */

/*
 * The PRP functions presented here are based on the paper:
 * Grantham, Jon. Frobenius Pseudoprimes. Math. Comp. 70 (2001), 873-891.
 */

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "mpz_aprcl32.h"


/* ******************************************************************
 * mpz_prp: (also called a Fermat pseudoprime)
 * A "pseudoprime" to the base a is a composite number n such that,
 * (a,n)=1 and a^(n-1) = 1 mod n
 * ******************************************************************/
int mpz_prp(mpz_t n, mpz_t a)
{
  mpz_t res;
  mpz_t nm1;

  if (mpz_cmp_ui(a, 2) < 0)
    return PRP_ERROR;

  if (mpz_cmp_ui(n, 2) < 0)
    return PRP_COMPOSITE;

  if (mpz_divisible_ui_p(n, 2))
  {
    if (mpz_cmp_ui(n, 2) == 0)
      return PRP_PRIME;
    else
      return PRP_COMPOSITE;
  }

  mpz_init_set_ui(res, 0);
  mpz_gcd(res, n, a);

  if (mpz_cmp_ui(res, 1) > 0)
  {
    mpz_clear(res);
    return PRP_COMPOSITE;
  }

  mpz_init_set(nm1, n);
  mpz_sub_ui(nm1, nm1, 1);
  mpz_powm(res, a, nm1, n);

  if (mpz_cmp_ui(res, 1) == 0)
  {
    mpz_clear(res);
    mpz_clear(nm1);
    return PRP_PRP;
  }
  else
  {
    mpz_clear(res);
    mpz_clear(nm1);
    return PRP_COMPOSITE;
  }

}/* method mpz_prp */


/* *************************************************************************
 * mpz_euler_prp: (also called a Solovay-Strassen pseudoprime)
 * An "Euler pseudoprime" to the base a is an odd composite number n with,
 * (a,n)=1 such that a^((n-1)/2)=(a/n) mod n [(a/n) is the Jacobi symbol]
 * *************************************************************************/
int mpz_euler_prp(mpz_t n, mpz_t a)
{
  mpz_t res;
  mpz_t exp;
  int ret = 0;

  if (mpz_cmp_ui(a, 2) < 0)
    return PRP_ERROR;

  if (mpz_cmp_ui(n, 2) < 0)
    return PRP_COMPOSITE;

  if (mpz_divisible_ui_p(n, 2))
  {
    if (mpz_cmp_ui(n, 2) == 0)
      return PRP_PRIME;
    else
      return PRP_COMPOSITE;
  }

  mpz_init_set_ui(res, 0);
  mpz_gcd(res, n, a);

  if (mpz_cmp_ui(res, 1) > 0)
  {
    mpz_clear(res);
    return PRP_COMPOSITE;
  }

  mpz_init_set(exp, n);

  mpz_sub_ui(exp, exp, 1); /* exp = n-1 */
  mpz_divexact_ui(exp, exp, 2); /* exp = (n-1)/2 */
  mpz_powm(res, a, exp, n);

  /* reuse exp to calculate jacobi(a,n) mod n */
  ret = mpz_jacobi(a,n);
  mpz_set(exp, n);
  if (ret == -1)
    mpz_sub_ui(exp, exp, 1);
  else if (ret == 1)
    mpz_add_ui(exp, exp, 1);
  mpz_mod(exp, exp, n);

  if (mpz_cmp(res, exp) == 0)
  {
    mpz_clear(res);
    mpz_clear(exp);
    return PRP_PRP;
  }
  else
  {
    mpz_clear(res);
    mpz_clear(exp);
    return PRP_COMPOSITE;
  }

}/* method mpz_euler_prp */


/* *********************************************************************************************
 * mpz_sprp: (also called a Miller-Rabin pseudoprime)
 * A "strong pseudoprime" to the base a is an odd composite n = (2^r)*s+1 with s odd such that
 * either a^s == 1 mod n, or a^((2^t)*s) == -1 mod n, for some integer t, with 0 <= t < r.
 * *********************************************************************************************/
int mpz_sprp(mpz_t n, mpz_t a)
{
  mpz_t s;
  mpz_t nm1;
  mpz_t mpz_test;
  unsigned long int r = 0;

  if (mpz_cmp_ui(a, 2) < 0)
    return PRP_ERROR;

  if (mpz_cmp_ui(n, 2) < 0)
    return PRP_COMPOSITE;

  if (mpz_divisible_ui_p(n, 2))
  {
    if (mpz_cmp_ui(n, 2) == 0)
      return PRP_PRIME;
    else
      return PRP_COMPOSITE;
  }

  mpz_init_set_ui(mpz_test, 0);
  mpz_init_set_ui(s, 0);
  mpz_init_set(nm1, n);
  mpz_sub_ui(nm1, nm1, 1);

  /***********************************************/
  /* Find s and r satisfying: n-1=(2^r)*s, s odd */
  r = mpz_scan1(nm1, 0);
  mpz_fdiv_q_2exp(s, nm1, r);


  /******************************************/
  /* Check a^((2^t)*s) mod n for 0 <= t < r */
  mpz_powm(mpz_test, a, s, n);
  if ( (mpz_cmp_ui(mpz_test, 1) == 0) || (mpz_cmp(mpz_test, nm1) == 0) )
  {
    mpz_clear(s);
    mpz_clear(nm1);
    mpz_clear(mpz_test);
    return PRP_PRP;
  }

  while ( --r )
  {
    /* mpz_test = mpz_test^2%n */
    mpz_mul(mpz_test, mpz_test, mpz_test);
    mpz_mod(mpz_test, mpz_test, n);

    if (mpz_cmp(mpz_test, nm1) == 0)
    {
      mpz_clear(s);
      mpz_clear(nm1);
      mpz_clear(mpz_test);
      return PRP_PRP;
    }
  }

  mpz_clear(s);
  mpz_clear(nm1);
  mpz_clear(mpz_test);
  return PRP_COMPOSITE;

}/* method mpz_sprp */


/* *************************************************************************
 * mpz_fibonacci_prp:
 * A "Fibonacci pseudoprime" with parameters (P,Q), P > 0, Q=+/-1, is a
 * composite n for which V_n == P mod n
 * [V is the Lucas V sequence with parameters P,Q]
 * *************************************************************************/
int mpz_fibonacci_prp(mpz_t n, long int p, long int q)
{
  mpz_t pmodn, zP;
  mpz_t vl, vh, ql, qh, tmp; /* used for calculating the Lucas V sequence */
  int s = 0, j = 0;

  if (p*p-4*q == 0)
    return PRP_ERROR;

  if (((q != 1) && (q != -1)) || (p <= 0))
    return PRP_ERROR;

  if (mpz_cmp_ui(n, 2) < 0)
    return PRP_COMPOSITE;

  if (mpz_divisible_ui_p(n, 2))
  {
    if (mpz_cmp_ui(n, 2) == 0)
      return PRP_PRIME;
    else
      return PRP_COMPOSITE;
  }

  mpz_init_set_ui(zP, p);
  mpz_init(pmodn);
  mpz_mod(pmodn, zP, n);

  /* mpz_lucasvmod(res, p, q, n, n); */
  mpz_init_set_si(vl, 2);
  mpz_init_set_si(vh, p);
  mpz_init_set_si(ql, 1);
  mpz_init_set_si(qh, 1);
  mpz_init_set_si(tmp,0);

  s = mpz_scan1(n, 0);
  for (j = mpz_sizeinbase(n,2)-1; j >= s+1; j--)
  {
    /* ql = ql*qh (mod n) */
    mpz_mul(ql, ql, qh);
    mpz_mod(ql, ql, n);
    if (mpz_tstbit(n,j) == 1)
    {
      /* qh = ql*q */
      mpz_mul_si(qh, ql, q);

      /* vl = vh*vl - p*ql (mod n) */
      mpz_mul(vl, vh, vl);
      mpz_mul_si(tmp, ql, p);
      mpz_sub(vl, vl, tmp);
      mpz_mod(vl, vl, n);

      /* vh = vh*vh - 2*qh (mod n) */
      mpz_mul(vh, vh, vh);
      mpz_mul_si(tmp, qh, 2);
      mpz_sub(vh, vh, tmp);
      mpz_mod(vh, vh, n);
    }
    else
    {
      /* qh = ql */
      mpz_set(qh, ql);

      /* vh = vh*vl - p*ql (mod n) */
      mpz_mul(vh, vh, vl);
      mpz_mul_si(tmp, ql, p);
      mpz_sub(vh, vh, tmp);
      mpz_mod(vh, vh, n);

      /* vl = vl*vl - 2*ql (mod n) */
      mpz_mul(vl, vl, vl);
      mpz_mul_si(tmp, ql, 2);
      mpz_sub(vl, vl, tmp);
      mpz_mod(vl, vl, n);
    }
  }
  /* ql = ql*qh */
  mpz_mul(ql, ql, qh);

  /* qh = ql*q */
  mpz_mul_si(qh, ql, q);

  /* vl = vh*vl - p*ql */
  mpz_mul(vl, vh, vl);
  mpz_mul_si(tmp, ql, p);
  mpz_sub(vl, vl, tmp);

  /* ql = ql*qh */
  mpz_mul(ql, ql, qh);

  for (j = 1; j <= s; j++)
  {
    /* vl = vl*vl - 2*ql (mod n) */
    mpz_mul(vl, vl, vl);
    mpz_mul_si(tmp, ql, 2);
    mpz_sub(vl, vl, tmp);
    mpz_mod(vl, vl, n);

    /* ql = ql*ql (mod n) */
    mpz_mul(ql, ql, ql);
    mpz_mod(ql, ql, n);
  }

  mpz_mod(vl, vl, n); /* vl contains our return value */

  if (mpz_cmp(vl, pmodn) == 0)
  {
    mpz_clear(zP);
    mpz_clear(pmodn);
    mpz_clear(vl);
    mpz_clear(vh);
    mpz_clear(ql);
    mpz_clear(qh);
    mpz_clear(tmp);
    return PRP_PRP;
  }
  mpz_clear(zP);
  mpz_clear(pmodn);
  mpz_clear(vl);
  mpz_clear(vh);
  mpz_clear(ql);
  mpz_clear(qh);
  mpz_clear(tmp);
  return PRP_COMPOSITE;

}/* method mpz_fibonacci_prp */


/* *******************************************************************************
 * mpz_lucas_prp:
 * A "Lucas pseudoprime" with parameters (P,Q) is a composite n with D=P^2-4Q,
 * (n,2QD)=1 such that U_(n-(D/n)) == 0 mod n [(D/n) is the Jacobi symbol]
 * *******************************************************************************/
int mpz_lucas_prp(mpz_t n, long int p, long int q)
{
  mpz_t zD;
  mpz_t res;
  mpz_t index;
  mpz_t uh, vl, vh, ql, qh, tmp; /* used for calculating the Lucas U sequence */
  int s = 0, j = 0;
  int ret = 0;
  long int d = p*p - 4*q;

  if (d == 0) /* Does not produce a proper Lucas sequence */
    return PRP_ERROR;

  if (mpz_cmp_ui(n, 2) < 0)
    return PRP_COMPOSITE;

  if (mpz_divisible_ui_p(n, 2))
  {
    if (mpz_cmp_ui(n, 2) == 0)
      return PRP_PRIME;
    else
      return PRP_COMPOSITE;
  }

  mpz_init(index);
  mpz_init_set_si(zD, d);
  mpz_init(res);

  mpz_mul_si(res, zD, q);
  mpz_mul_ui(res, res, 2);
  mpz_gcd(res, res, n);
  if ((mpz_cmp(res, n) != 0) && (mpz_cmp_ui(res, 1) > 0))
  {
    mpz_clear(zD);
    mpz_clear(res);
    mpz_clear(index);
    return PRP_COMPOSITE;
  }

  /* index = n-(D/n), where (D/n) is the Jacobi symbol */
  mpz_set(index, n);
  ret = mpz_jacobi(zD, n);
  if (ret == -1)
    mpz_add_ui(index, index, 1);
  else if (ret == 1)
    mpz_sub_ui(index, index, 1);

  /* mpz_lucasumod(res, p, q, index, n); */
  mpz_init_set_si(uh, 1);
  mpz_init_set_si(vl, 2);
  mpz_init_set_si(vh, p);
  mpz_init_set_si(ql, 1);
  mpz_init_set_si(qh, 1);
  mpz_init_set_si(tmp,0);

  s = mpz_scan1(index, 0);
  for (j = mpz_sizeinbase(index,2)-1; j >= s+1; j--)
  {
    /* ql = ql*qh (mod n) */
    mpz_mul(ql, ql, qh);
    mpz_mod(ql, ql, n);
    if (mpz_tstbit(index,j) == 1)
    {
      /* qh = ql*q */
      mpz_mul_si(qh, ql, q);

      /* uh = uh*vh (mod n) */
      mpz_mul(uh, uh, vh);
      mpz_mod(uh, uh, n);

      /* vl = vh*vl - p*ql (mod n) */
      mpz_mul(vl, vh, vl);
      mpz_mul_si(tmp, ql, p);
      mpz_sub(vl, vl, tmp);
      mpz_mod(vl, vl, n);

      /* vh = vh*vh - 2*qh (mod n) */
      mpz_mul(vh, vh, vh);
      mpz_mul_si(tmp, qh, 2);
      mpz_sub(vh, vh, tmp);
      mpz_mod(vh, vh, n);
    }
    else
    {
      /* qh = ql */
      mpz_set(qh, ql);

      /* uh = uh*vl - ql (mod n) */
      mpz_mul(uh, uh, vl);
      mpz_sub(uh, uh, ql);
      mpz_mod(uh, uh, n);

      /* vh = vh*vl - p*ql (mod n) */
      mpz_mul(vh, vh, vl);
      mpz_mul_si(tmp, ql, p);
      mpz_sub(vh, vh, tmp);
      mpz_mod(vh, vh, n);

      /* vl = vl*vl - 2*ql (mod n) */
      mpz_mul(vl, vl, vl);
      mpz_mul_si(tmp, ql, 2);
      mpz_sub(vl, vl, tmp);
      mpz_mod(vl, vl, n);
    }
  }
  /* ql = ql*qh */
  mpz_mul(ql, ql, qh);

  /* qh = ql*q */
  mpz_mul_si(qh, ql, q);

  /* uh = uh*vl - ql */
  mpz_mul(uh, uh, vl);
  mpz_sub(uh, uh, ql);

  /* vl = vh*vl - p*ql */
  mpz_mul(vl, vh, vl);
  mpz_mul_si(tmp, ql, p);
  mpz_sub(vl, vl, tmp);

  /* ql = ql*qh */
  mpz_mul(ql, ql, qh);

  for (j = 1; j <= s; j++)
  {
    /* uh = uh*vl (mod n) */
    mpz_mul(uh, uh, vl);
    mpz_mod(uh, uh, n);

    /* vl = vl*vl - 2*ql (mod n) */
    mpz_mul(vl, vl, vl);
    mpz_mul_si(tmp, ql, 2);
    mpz_sub(vl, vl, tmp);
    mpz_mod(vl, vl, n);

    /* ql = ql*ql (mod n) */
    mpz_mul(ql, ql, ql);
    mpz_mod(ql, ql, n);
  }

  mpz_mod(res, uh, n); /* uh contains our return value */

  mpz_clear(zD);
  mpz_clear(index);
  mpz_clear(uh);
  mpz_clear(vl);
  mpz_clear(vh);
  mpz_clear(ql);
  mpz_clear(qh);
  mpz_clear(tmp);

  if (mpz_cmp_ui(res, 0) == 0)
  {
    mpz_clear(res);
    return PRP_PRP;
  }
  else
  {
    mpz_clear(res);
    return PRP_COMPOSITE;
  }

}/* method mpz_lucas_prp */


/* *********************************************************************************************
 * mpz_stronglucas_prp:
 * A "strong Lucas pseudoprime" with parameters (P,Q) is a composite n = (2^r)*s+(D/n), where
 * s is odd, D=P^2-4Q, and (n,2QD)=1 such that either U_s == 0 mod n or V_((2^t)*s) == 0 mod n
 * for some t, 0 <= t < r. [(D/n) is the Jacobi symbol]
 * *********************************************************************************************/
int mpz_stronglucas_prp(mpz_t n, long int p, long int q)
{
  mpz_t zD;
  mpz_t s;
  mpz_t nmj; /* n minus jacobi(D/n) */
  mpz_t res;
  mpz_t uh, vl, vh, ql, qh, tmp; /* these are needed for the LucasU and LucasV part of this function */
  long int d = p*p - 4*q;
  unsigned long int r = 0;
  int ret = 0;
  unsigned long int j = 0;

  if (d == 0) /* Does not produce a proper Lucas sequence */
    return PRP_ERROR;

  if (mpz_cmp_ui(n, 2) < 0)
    return PRP_COMPOSITE;

  if (mpz_divisible_ui_p(n, 2))
  {
    if (mpz_cmp_ui(n, 2) == 0)
      return PRP_PRIME;
    else
      return PRP_COMPOSITE;
  }

  mpz_init_set_si(zD, d);
  mpz_init(res);

  mpz_mul_si(res, zD, q);
  mpz_mul_ui(res, res, 2);
  mpz_gcd(res, res, n);
  if ((mpz_cmp(res, n) != 0) && (mpz_cmp_ui(res, 1) > 0))
  {
    mpz_clear(zD);
    mpz_clear(res);
    return PRP_COMPOSITE;
  }

  mpz_init(s);
  mpz_init(nmj);

  /* nmj = n - (D/n), where (D/n) is the Jacobi symbol */
  mpz_set(nmj, n);
  ret = mpz_jacobi(zD, n);
  if (ret == -1)
    mpz_add_ui(nmj, nmj, 1);
  else if (ret == 1)
    mpz_sub_ui(nmj, nmj, 1);

  r = mpz_scan1(nmj, 0);
  mpz_fdiv_q_2exp(s, nmj, r);

  /* make sure U_s == 0 mod n or V_((2^t)*s) == 0 mod n, for some t, 0 <= t < r */
  mpz_init_set_si(uh, 1);
  mpz_init_set_si(vl, 2);
  mpz_init_set_si(vh, p);
  mpz_init_set_si(ql, 1);
  mpz_init_set_si(qh, 1);
  mpz_init_set_si(tmp,0);

  for (j = mpz_sizeinbase(s,2)-1; j >= 1; j--)
  {
    /* ql = ql*qh (mod n) */
    mpz_mul(ql, ql, qh);
    mpz_mod(ql, ql, n);
    if (mpz_tstbit(s,j) == 1)
    {
      /* qh = ql*q */
      mpz_mul_si(qh, ql, q);

      /* uh = uh*vh (mod n) */
      mpz_mul(uh, uh, vh);
      mpz_mod(uh, uh, n);

      /* vl = vh*vl - p*ql (mod n) */
      mpz_mul(vl, vh, vl);
      mpz_mul_si(tmp, ql, p);
      mpz_sub(vl, vl, tmp);
      mpz_mod(vl, vl, n);

      /* vh = vh*vh - 2*qh (mod n) */
      mpz_mul(vh, vh, vh);
      mpz_mul_si(tmp, qh, 2);
      mpz_sub(vh, vh, tmp);
      mpz_mod(vh, vh, n);
    }
    else
    {
      /* qh = ql */
      mpz_set(qh, ql);

      /* uh = uh*vl - ql (mod n) */
      mpz_mul(uh, uh, vl);
      mpz_sub(uh, uh, ql);
      mpz_mod(uh, uh, n);

      /* vh = vh*vl - p*ql (mod n) */
      mpz_mul(vh, vh, vl);
      mpz_mul_si(tmp, ql, p);
      mpz_sub(vh, vh, tmp);
      mpz_mod(vh, vh, n);

      /* vl = vl*vl - 2*ql (mod n) */
      mpz_mul(vl, vl, vl);
      mpz_mul_si(tmp, ql, 2);
      mpz_sub(vl, vl, tmp);
      mpz_mod(vl, vl, n);
    }
  }
  /* ql = ql*qh */
  mpz_mul(ql, ql, qh);

  /* qh = ql*q */
  mpz_mul_si(qh, ql, q);

  /* uh = uh*vl - ql */
  mpz_mul(uh, uh, vl);
  mpz_sub(uh, uh, ql);

  /* vl = vh*vl - p*ql */
  mpz_mul(vl, vh, vl);
  mpz_mul_si(tmp, ql, p);
  mpz_sub(vl, vl, tmp);

  /* ql = ql*qh */
  mpz_mul(ql, ql, qh);

  mpz_mod(uh, uh, n);
  mpz_mod(vl, vl, n);

  /* uh contains LucasU_s and vl contains LucasV_s */
  if ((mpz_cmp_ui(uh, 0) == 0) || (mpz_cmp_ui(vl, 0) == 0))
  {
    mpz_clear(zD);
    mpz_clear(s);
    mpz_clear(nmj);
    mpz_clear(res);
    mpz_clear(uh);
    mpz_clear(vl);
    mpz_clear(vh);
    mpz_clear(ql);
    mpz_clear(qh);
    mpz_clear(tmp);
    return PRP_PRP;
  }

  for (j = 1; j < r; j++)
  {
    /* vl = vl*vl - 2*ql (mod n) */
    mpz_mul(vl, vl, vl);
    mpz_mul_si(tmp, ql, 2);
    mpz_sub(vl, vl, tmp);
    mpz_mod(vl, vl, n);

    /* ql = ql*ql (mod n) */
    mpz_mul(ql, ql, ql);
    mpz_mod(ql, ql, n);

    if (mpz_cmp_ui(vl, 0) == 0)
    {
      mpz_clear(zD);
      mpz_clear(s);
      mpz_clear(nmj);
      mpz_clear(res);
      mpz_clear(uh);
      mpz_clear(vl);
      mpz_clear(vh);
      mpz_clear(ql);
      mpz_clear(qh);
      mpz_clear(tmp);
      return PRP_PRP;
    }
  }

  mpz_clear(zD);
  mpz_clear(s);
  mpz_clear(nmj);
  mpz_clear(res);
  mpz_clear(uh);
  mpz_clear(vl);
  mpz_clear(vh);
  mpz_clear(ql);
  mpz_clear(qh);
  mpz_clear(tmp);
  return PRP_COMPOSITE;

}/* method mpz_stronglucas_prp */


/* *******************************************************************************************
 * mpz_extrastronglucas_prp:
 * Let U_n = LucasU(p,1), V_n = LucasV(p,1), and D=p^2-4.
 * An "extra strong Lucas pseudoprime" to the base p is a composite n = (2^r)*s+(D/n), where
 * s is odd and (n,2D)=1, such that either U_s == 0 mod n and V_s == +/-2 mod n, or
 * V_((2^t)*s) == 0 mod n for some t with 0 <= t < r-1 [(D/n) is the Jacobi symbol]
 * *******************************************************************************************/
int mpz_extrastronglucas_prp(mpz_t n, long int p)
{
  mpz_t zD;
  mpz_t s;
  mpz_t nmj; /* n minus jacobi(D/n) */
  mpz_t res;
  mpz_t uh, vl, vh, ql, qh, tmp; /* these are needed for the LucasU and LucasV part of this function */
  long int d = p*p - 4;
  long int q = 1;
  unsigned long int r = 0;
  int ret = 0;
  unsigned long int j = 0;

  if (d == 0) /* Does not produce a proper Lucas sequence */
    return PRP_ERROR;

  if (mpz_cmp_ui(n, 2) < 0)
    return PRP_COMPOSITE;

  if (mpz_divisible_ui_p(n, 2))
  {
    if (mpz_cmp_ui(n, 2) == 0)
      return PRP_PRIME;
    else
      return PRP_COMPOSITE;
  }

  mpz_init_set_si(zD, d);
  mpz_init(res);

  mpz_mul_ui(res, zD, 2);
  mpz_gcd(res, res, n);
  if ((mpz_cmp(res, n) != 0) && (mpz_cmp_ui(res, 1) > 0))
  {
    mpz_clear(zD);
    mpz_clear(res);
    return PRP_COMPOSITE;
  }

  mpz_init(s);
  mpz_init(nmj);

  /* nmj = n - (D/n), where (D/n) is the Jacobi symbol */
  mpz_set(nmj, n);
  ret = mpz_jacobi(zD, n);
  if (ret == -1)
    mpz_add_ui(nmj, nmj, 1);
  else if (ret == 1)
    mpz_sub_ui(nmj, nmj, 1);

  r = mpz_scan1(nmj, 0);
  mpz_fdiv_q_2exp(s, nmj, r);

  /* make sure that either (U_s == 0 mod n and V_s == +/-2 mod n), or */
  /* V_((2^t)*s) == 0 mod n for some t with 0 <= t < r-1           */
  mpz_init_set_si(uh, 1);
  mpz_init_set_si(vl, 2);
  mpz_init_set_si(vh, p);
  mpz_init_set_si(ql, 1);
  mpz_init_set_si(qh, 1);
  mpz_init_set_si(tmp,0);

  for (j = mpz_sizeinbase(s,2)-1; j >= 1; j--)
  {
    /* ql = ql*qh (mod n) */
    mpz_mul(ql, ql, qh);
    mpz_mod(ql, ql, n);
    if (mpz_tstbit(s,j) == 1)
    {
      /* qh = ql*q */
      mpz_mul_si(qh, ql, q);

      /* uh = uh*vh (mod n) */
      mpz_mul(uh, uh, vh);
      mpz_mod(uh, uh, n);

      /* vl = vh*vl - p*ql (mod n) */
      mpz_mul(vl, vh, vl);
      mpz_mul_si(tmp, ql, p);
      mpz_sub(vl, vl, tmp);
      mpz_mod(vl, vl, n);

      /* vh = vh*vh - 2*qh (mod n) */
      mpz_mul(vh, vh, vh);
      mpz_mul_si(tmp, qh, 2);
      mpz_sub(vh, vh, tmp);
      mpz_mod(vh, vh, n);
    }
    else
    {
      /* qh = ql */
      mpz_set(qh, ql);

      /* uh = uh*vl - ql (mod n) */
      mpz_mul(uh, uh, vl);
      mpz_sub(uh, uh, ql);
      mpz_mod(uh, uh, n);

      /* vh = vh*vl - p*ql (mod n) */
      mpz_mul(vh, vh, vl);
      mpz_mul_si(tmp, ql, p);
      mpz_sub(vh, vh, tmp);
      mpz_mod(vh, vh, n);

      /* vl = vl*vl - 2*ql (mod n) */
      mpz_mul(vl, vl, vl);
      mpz_mul_si(tmp, ql, 2);
      mpz_sub(vl, vl, tmp);
      mpz_mod(vl, vl, n);
    }
  }
  /* ql = ql*qh */
  mpz_mul(ql, ql, qh);

  /* qh = ql*q */
  mpz_mul_si(qh, ql, q);

  /* uh = uh*vl - ql */
  mpz_mul(uh, uh, vl);
  mpz_sub(uh, uh, ql);

  /* vl = vh*vl - p*ql */
  mpz_mul(vl, vh, vl);
  mpz_mul_si(tmp, ql, p);
  mpz_sub(vl, vl, tmp);

  /* ql = ql*qh */
  mpz_mul(ql, ql, qh);

  mpz_mod(uh, uh, n);
  mpz_mod(vl, vl, n);

  /* tmp = n-2, for the following comparison */
  mpz_sub_ui(tmp, n, 2);

  /* uh contains LucasU_s and vl contains LucasV_s */
  if ((mpz_cmp_ui(uh, 0) == 0) && ((mpz_cmp(vl, tmp) == 0) || (mpz_cmp_si(vl, 2) == 0)))
  {
    mpz_clear(zD);
    mpz_clear(s);
    mpz_clear(nmj);
    mpz_clear(res);
    mpz_clear(uh);
    mpz_clear(vl);
    mpz_clear(vh);
    mpz_clear(ql);
    mpz_clear(qh);
    mpz_clear(tmp);
    return PRP_PRP;
  }

  for (j = 1; j < r-1; j++)
  {
    /* vl = vl*vl - 2*ql (mod n) */
    mpz_mul(vl, vl, vl);
    mpz_mul_si(tmp, ql, 2);
    mpz_sub(vl, vl, tmp);
    mpz_mod(vl, vl, n);

    /* ql = ql*ql (mod n) */
    mpz_mul(ql, ql, ql);
    mpz_mod(ql, ql, n);

    if (mpz_cmp_ui(vl, 0) == 0)
    {
      mpz_clear(zD);
      mpz_clear(s);
      mpz_clear(nmj);
      mpz_clear(res);
      mpz_clear(uh);
      mpz_clear(vl);
      mpz_clear(vh);
      mpz_clear(ql);
      mpz_clear(qh);
      mpz_clear(tmp);
      return PRP_PRP;
    }
  }

  mpz_clear(zD);
  mpz_clear(s);
  mpz_clear(nmj);
  mpz_clear(res);
  mpz_clear(uh);
  mpz_clear(vl);
  mpz_clear(vh);
  mpz_clear(ql);
  mpz_clear(qh);
  mpz_clear(tmp);
  return PRP_COMPOSITE;

}/* method mpz_extrastronglucas_prp */


/* ***********************************************************************************************
 * mpz_selfridge_prp:
 * A "Lucas-Selfridge pseudoprime" n is a "Lucas pseudoprime" using Selfridge parameters of:
 * Find the first element D in the sequence {5, -7, 9, -11, 13, ...} such that Jacobi(D,n) = -1
 * Then use P=1 and Q=(1-D)/4 in the Lucas pseudoprime test.
 * Make sure n is not a perfect square, otherwise the search for D will only stop when D=n.
 * ***********************************************************************************************/
int mpz_selfridge_prp(mpz_t n)
{
  long int d = 5, p = 1, q = 0;
  int max_d = 1000000;
  int jacobi = 0;
  mpz_t zD;

  if (mpz_cmp_ui(n, 2) < 0)
    return PRP_COMPOSITE;

  if (mpz_divisible_ui_p(n, 2))
  {
    if (mpz_cmp_ui(n, 2) == 0)
      return PRP_PRIME;
    else
      return PRP_COMPOSITE;
  }

  mpz_init_set_ui(zD, d);

  while (1)
  {
    jacobi = mpz_jacobi(zD, n);

    /* if jacobi == 0, d is a factor of n, therefore n is composite... */
    /* if d == n, then either n is either prime or 9... */
    if (jacobi == 0)
    {
      if ((mpz_cmpabs(zD, n) == 0) && (mpz_cmp_ui(zD, 9) != 0))
      {
        mpz_clear(zD);
        return PRP_PRIME;
      }
      else
      {
        mpz_clear(zD);
        return PRP_COMPOSITE;
      }
    }
    if (jacobi == -1)
      break;

    /* if we get to the 5th d, make sure we aren't dealing with a square... */
    if (d == 13)
    {
      if (mpz_perfect_square_p(n))
      {
        mpz_clear(zD);
        return PRP_COMPOSITE;
      }
    }

    if (d < 0)
    {
      d *= -1;
      d += 2;
    }
    else
    {
      d += 2;
      d *= -1;
    }

    /* make sure we don't search forever */
    if (d >= max_d)
    {
      mpz_clear(zD);
      return PRP_ERROR;
    }

    mpz_set_si(zD, d);
  }
  mpz_clear(zD);

  q = (1-d)/4;

  return mpz_lucas_prp(n, p, q);

}/* method mpz_selfridge_prp */


/* *********************************************************************************************************
 * mpz_strongselfridge_prp:
 * A "strong Lucas-Selfridge pseudoprime" n is a "strong Lucas pseudoprime" using Selfridge parameters of:
 * Find the first element D in the sequence {5, -7, 9, -11, 13, ...} such that Jacobi(D,n) = -1
 * Then use P=1 and Q=(1-D)/4 in the strong Lucase pseudoprime test.
 * Make sure n is not a perfect square, otherwise the search for D will only stop when D=n.
 * **********************************************************************************************************/
int mpz_strongselfridge_prp(mpz_t n)
{
  long int d = 5, p = 1, q = 0;
  int max_d = 1000000;
  int jacobi = 0;
  mpz_t zD;

  if (mpz_cmp_ui(n, 2) < 0)
    return PRP_COMPOSITE;

  if (mpz_divisible_ui_p(n, 2))
  {
    if (mpz_cmp_ui(n, 2) == 0)
      return PRP_PRIME;
    else
      return PRP_COMPOSITE;
  }

  mpz_init_set_ui(zD, d);

  while (1)
  {
    jacobi = mpz_jacobi(zD, n);

    /* if jacobi == 0, d is a factor of n, therefore n is composite... */
    /* if d == n, then either n is either prime or 9... */
    if (jacobi == 0)
    {
      if ((mpz_cmpabs(zD, n) == 0) && (mpz_cmp_ui(zD, 9) != 0))
      {
        mpz_clear(zD);
        return PRP_PRIME;
      }
      else
      {
        mpz_clear(zD);
        return PRP_COMPOSITE;
      }
    }
    if (jacobi == -1)
      break;

    /* if we get to the 5th d, make sure we aren't dealing with a square... */
    if (d == 13)
    {
      if (mpz_perfect_square_p(n))
      {
        mpz_clear(zD);
        return PRP_COMPOSITE;
      }
    }

    if (d < 0)
    {
      d *= -1;
      d += 2;
    }
    else
    {
      d += 2;
      d *= -1;
    }

    /* make sure we don't search forever */
    if (d >= max_d)
    {
      mpz_clear(zD);
      return PRP_ERROR;
    }

    mpz_set_si(zD, d);
  }
  mpz_clear(zD);

  q = (1-d)/4;

  return mpz_stronglucas_prp(n, p, q);

}/* method mpz_strongselfridge_prp */


/* **********************************************************************************
 * mpz_bpsw_prp:
 * A "Baillie-Pomerance-Selfridge-Wagstaff pseudoprime" is a composite n such that
 * n is a strong pseudoprime to the base 2 and
 * n is a Lucas pseudoprime using the Selfridge parameters.
 * **********************************************************************************/
int mpz_bpsw_prp(mpz_t n)
{
  int ret = 0;
  mpz_t two;

  mpz_init_set_ui(two, 2);

  ret = mpz_sprp(n, two);
  mpz_clear(two);

  /* with a base of 2, mpz_sprp, won't return PRP_ERROR */
  /* so, only check for PRP_COMPOSITE or PRP_PRIME here */
  if ((ret == PRP_COMPOSITE) || (ret == PRP_PRIME))
    return ret;

  return mpz_selfridge_prp(n);

}/* method mpz_bpsw_prp */


/* ****************************************************************************************
 * mpz_strongbpsw_prp:
 * A "strong Baillie-Pomerance-Selfridge-Wagstaff pseudoprime" is a composite n such that
 * n is a strong pseudoprime to the base 2 and
 * n is a strong Lucas pseudoprime using the Selfridge parameters.
 * ****************************************************************************************/
int mpz_strongbpsw_prp(mpz_t n)
{
  int ret = 0;
  mpz_t two;

  mpz_init_set_ui(two, 2);

  ret = mpz_sprp(n, two);
  mpz_clear(two);

  /* with a base of 2,  mpz_sprp won't return PRP_ERROR */
  /* so, only check for PRP_COMPOSITE or PRP_PRIME here */
  if ((ret == PRP_COMPOSITE) || (ret == PRP_PRIME))
    return ret;

  return mpz_strongselfridge_prp(n);

}/* method mpz_strongbpsw_prp */


/* ************************
 * mpz_frobenius_prp
 * ************************/

/* ************************
 * mpz_strongfrobenius_prp
 * ************************/

/* ************************
 * mpz_elliptic_prp
 * ************************/


/********************************************/
/********************************************/
/* End of probable primality (prp) routines */
/********************************************/
/********************************************/


/***********************************/
/***********************************/
/* Start of prime proving routines */
/***********************************/
/***********************************/

/* **********************************************************************************
 * APR-CL (also known as APRT-CLE) is a prime proving algorithm developed by:
 * L. Adleman, C. Pomerance, R. Rumely, H. Cohen, and H. W. Lenstra
 * APRT-CLE = Adleman-Pomerance-Rumely Test Cohen-Lenstra Extended version
 * You can find all the details of this implementation in the Cohen & Lenstra paper:
 *    H. Cohen and A. K. Lenstra, "Implementation of a new primality test",
 *    Math. Comp., 48 (1987) 103--121
 *
 * ----------------------------------------------------------------------------------
 *
 * This C/GMP version is a conversion of Dario Alpern's Java based APRT-CLE code
 * His code was based on Yuji Kida's UBASIC code
 *
 * Based on APRT-CLE Written by Dario Alejandro Alpern (Buenos Aires - Argentina)
 * From: Last updated September 10th, 2011. See http://www.alpertron.com.ar/ECM.HTM
 *
 * On 2012/11/12 Dario Alpern has approved the conversion, from Java to C/GMP, of 
 * his implementation of the APR-CL algorithm, and that it be licensed under the LGPL.
 *
 * ----------------------------------------------------------------------------------
 *
 * With improvements based on Jason Moxham's APRCL v1.15 code, from 2003/01/01
 *
 * On 2013/04/14 Toby Moxham has approved the APR-CL code and data tables, 
 * originally written by his brother Jason Moxham on 2003/01/01, to be released 
 * under the LGPL.
 *
 * *********************************************************************************/

/* int PWmax = 32, Qmax = 55441, LEVELmax = 11; */
/* QMax is the largest Q in the aiQ array */
/* PWmax is the largest p^k that divides any (Q-1) */
/* #define Qmax 55441 */
/* #define PWmax 32 */
/*
 2-max =  2^5 =  32 from t[ 4]=    4324320
 3-max =  3^3 =  27 from t[ 4]=    4324320
 5-max =  5^2 =  25 from t[ 6]=  367567200
 7-max =  7^1 =   7 from t[ 1]=       5040
11-max = 11^1 =  11 from t[ 2]=      55440
13-max = 13^1 =  13 from t[ 3]=     720720
17-max = 17^1 =  17 from t[ 5]=   73513440
19-max = 19^1 =  19 from t[ 7]= 1396755360
*/

/* largest value of p^k that divides any t value */
#define PWmax 32

/* Qmax[]: list of largest q-prime for each t value */
/* this determines how much memory the program will need at run-time */
/* memory requirements will be around Qmax*4 bytes */
int Qmax[] = {61,2521,55441,180181,4324321,10501921,367567201,232792561};

/* number of values in Qmax[], aiNP[], aiNQ[] */
int LEVELmax = 8;

/* list of primes that divide our t values */
int aiP[] = {2,3,5,7,11,13,17,19};

/* list of q-primes from all t values */
int aiQ[] = {2,3,5,7,11,13,31,61,17,19,29,37,41,43,71,73,113,127,181,211,241,
281,337,421,631,1009,2521,23,67,89,199,331,397,463,617,661,881,991,1321,2311,
3697,4621,9241,18481,55441,53,79,131,157,313,521,547,859,911,937,1093,1171,
1873,2003,2341,2731,2861,3121,3433,6007,6553,8009,8191,8581,16381,20021,20593,
21841,25741,36037,48049,51481,65521,72073,120121,180181,97,109,271,353,379,433,
541,673,757,1249,2017,2081,2161,2377,2971,3169,3361,3511,4159,5281,7393,7561,
7723,8317,8737,9829,13729,14561,15121,16633,23761,24571,26209,28081,30241,
38611,39313,47521,66529,96097,108109,110881,123553,131041,196561,216217,270271,
332641,393121,432433,540541,617761,4324321,103,137,239,307,409,443,613,919,953,
1021,1123,1327,1361,1429,1531,1871,2143,2381,2857,3061,3571,3673,4421,4591,
5237,6121,6427,6733,7481,8161,9181,9283,9521,10099,10711,12241,12377,12853,
14281,15913,16831,17137,17681,19891,22441,23563,23869,24481,27847,29173,29921,
30941,34273,36721,42841,43759,46411,47737,52361,53857,59671,63649,70687,72931,
74257,78541,79561,87517,92821,97241,100981,102103,116689,117811,128521,145861,
148513,157081,161569,167077,185641,201961,209441,235621,238681,269281,291721,
314161,371281,388961,417691,445537,471241,477361,514081,565489,612613,656371,
680681,700129,816817,1633633,1670761,1837837,2625481,4084081,5250961,5654881,
8168161,9189181,10501921,101,151,401,601,701,1051,1201,1301,1801,1951,2551,
2801,3301,3851,4201,4951,5101,5851,6301,7151,9901,11551,11701,12601,14851,
15401,15601,17551,17851,18701,19801,21601,23801,28051,33151,34651,40801,42901,
44201,50051,53551,54601,56101,66301,70201,77351,79201,81901,91801,92401,93601,
103951,107101,109201,118801,122401,140401,150151,151201,160651,193051,198901,
200201,218401,224401,232051,243101,257401,300301,321301,367201,415801,428401,
448801,450451,504901,530401,600601,673201,729301,795601,800801,982801,1029601,
1093951,1178101,1201201,1458601,2088451,2187901,2402401,2570401,2702701,
3088801,3141601,3712801,5105101,5834401,6806801,7068601,8353801,17503201,
22972951,52509601,183783601,367567201,191,229,419,457,571,647,761,1483,1597,
2053,2129,2281,2927,3041,3877,4447,4523,4561,4789,6271,6689,6841,6917,7411,
7753,8209,8779,8893,10337,11287,11971,12541,13339,13567,13681,14821,16417,
17291,17443,18089,19381,20521,20749,21319,21737,22573,25841,27361,28729,29641,
30097,31123,35531,35569,35911,38039,39521,40699,43891,46817,47881,48907,51871,
53353,56431,57457,58787,59281,63841,71821,72353,75583,77521,87211,90289,97813,
105337,106591,108529,114913,117041,124489,131671,134369,135661,139537,140449,
146719,163021,177841,186733,207481,213181,217057,217361,225721,251941,279073,
287281,300961,302329,342343,351121,377911,391249,406981,451441,456457,461891,
489061,511633,526681,554269,568481,608609,651169,652081,697681,733591,782497,
790021,813961,895357,1027027,1053361,1058149,1108537,1133731,1264033,1279081,
1369369,1492261,1580041,1790713,1813969,1867321,1939939,2217073,2238391,
2282281,2351441,2489761,2645371,2771341,2934361,2984521,3233231,3627937,
3837241,3912481,3979361,4157011,4232593,4476781,5135131,5372137,5868721,
6046561,6348889,6651217,6715171,6846841,7162849,7674481,9767521,11737441,
12471031,12697777,17907121,24942061,27387361,31744441,35814241,41081041,
46558513,53721361,107442721,174594421,232792561};
/* number of q primes in the above array: 530 */

/* a primitive root for each q-prime in the above array */
int aiG[] = {1,2,2,3,2,2,3,2,3,2,2,2,6,3,7,5,3,3,2,2,7,3,10,2,3,11,17,5,2,3,3,
3,5,3,3,2,3,6,13,3,5,2,13,13,38,2,3,2,5,10,3,2,2,17,5,5,2,10,5,7,3,2,7,5,3,10,
3,17,6,2,3,5,11,6,2,17,17,17,5,29,6,5,6,6,3,2,5,2,5,2,7,5,3,23,5,10,7,22,7,3,7,
5,13,3,6,5,10,23,6,11,15,7,7,11,19,11,3,10,17,19,5,2,69,5,17,11,15,7,19,23,5,
14,19,17,5,3,7,5,21,2,2,7,3,10,2,3,3,6,2,14,3,3,11,6,2,5,3,11,3,7,3,2,6,7,2,2,
3,2,3,7,6,5,19,5,6,5,3,2,14,2,2,11,6,2,3,2,5,37,23,3,3,5,6,5,12,7,5,10,5,2,7,2,
6,3,7,5,7,2,22,2,5,26,13,5,41,13,3,10,29,7,14,37,19,3,12,29,19,14,33,13,2,2,6,
28,5,5,19,5,7,19,29,13,23,6,7,2,6,3,7,2,7,11,2,11,3,6,3,6,2,11,6,6,2,10,7,2,7,
6,11,2,6,23,3,2,2,13,7,3,2,3,2,13,6,6,2,3,22,7,22,14,14,13,2,13,34,22,6,6,17,
17,7,11,6,17,2,2,2,3,22,31,3,2,29,2,2,11,31,19,13,2,2,23,37,23,7,19,3,17,7,3,6,
31,23,7,2,19,26,6,31,26,29,10,14,3,19,29,37,6,37,41,23,19,6,2,13,3,5,6,2,11,2,
3,7,5,3,2,3,5,11,2,11,3,22,2,2,10,7,11,5,3,3,10,14,2,3,22,2,10,6,2,3,7,11,2,14,
6,6,3,7,22,7,10,3,6,11,12,7,3,2,3,3,29,2,7,7,3,5,2,7,17,2,3,5,7,13,23,2,10,3,
23,5,3,11,3,3,10,5,17,6,6,7,5,31,10,10,6,17,6,10,13,7,7,3,29,3,7,6,29,5,18,17,
13,29,6,3,3,22,14,14,6,10,17,13,6,7,34,2,5,2,10,31,43,6,13,13,21,29,2,5,7,17,3,
22,7,7,7,29,14,5,13,21,6,10,15,6,2,5,14,14,11,5,7,23,13,7,37,29,11,5,13,22,37,
58,26,29,5,43,23,2,71};

/* number of primes from aiP, not necessarily in order, that divides each t */
int aiNP[] = {3,4,5,6,6,7,7,8};

/* number of q-primes for each t */
int aiNQ[] = {8,27,45,81,134,245,351,424};

/*         t     |       e(t)   | #Qp |      Qmax  |   divisors of t         */
/* --------------|--------------|-----|------------|------------------------ */
int aiT[] =  {
          60, /* | 6.8144 E   9 |   8 |         61 | p={2,3,5}               */
        5040, /* | 1.5321 E  52 |  27 |       2521 | p={2,3,5,7}             */
       55440, /* | 4.9209 E 106 |  45 |      55441 | p={2,3,5,7,11}          */
      720720, /* | 2.5992 E 237 |  81 |     180181 | p={2,3,5,7,11,13}       */
     4324320, /* | 7.9285 E 455 | 134 |    4324321 | p={2,3,5,7,11,13}       */
    73513440, /* | 7.0821 E 966 | 245 |   10501921 | p={2,3,5,7,11,13,17}    */
   367567200, /* | 6.2087 E1501 | 351 |  367567201 | p={2,3,5,7,11,13,17}    */
  1396755360};/* | 4.0165 E1913 | 424 |  232792561 | p={2,3,5,7,11,13,17,19} */


int aiInv[PWmax];
mpz_t biTmp;
mpz_t biExp;
mpz_t biN;
mpz_t biR;
mpz_t biS;
mpz_t biT;
mpz_t biU;
mpz_t *aiJS; /* [PWmax] */
mpz_t *aiJW; /* [PWmax] */
mpz_t *aiJX; /* [PWmax] */
mpz_t *aiJ0; /* [PWmax] */
mpz_t *aiJ1; /* [PWmax] */
mpz_t *aiJ2; /* [PWmax] */
mpz_t *aiJ00; /* [PWmax] */
mpz_t *aiJ01; /* [PWmax] */
int NumberLength; /* Length of multiple precision nbrs */
mpz_t TestNbr;

/* ============================================================================================== */

void allocate_vars()
{
  int i = 0;
  aiJS = malloc(PWmax * sizeof(mpz_t));
  aiJW = malloc(PWmax * sizeof(mpz_t));
  aiJX = malloc(PWmax * sizeof(mpz_t));
  aiJ0 = malloc(PWmax * sizeof(mpz_t));
  aiJ1 = malloc(PWmax * sizeof(mpz_t));
  aiJ2 = malloc(PWmax * sizeof(mpz_t));
  aiJ00 = malloc(PWmax * sizeof(mpz_t));
  aiJ01 = malloc(PWmax * sizeof(mpz_t));
  for (i = 0 ; i < PWmax; i++)
  {
    mpz_init(aiJS[i]);
    mpz_init(aiJW[i]);
    mpz_init(aiJX[i]);
    mpz_init(aiJ0[i]);
    mpz_init(aiJ1[i]);
    mpz_init(aiJ2[i]);
    mpz_init(aiJ00[i]);
    mpz_init(aiJ01[i]);
  }

  mpz_init(TestNbr);
  mpz_init(biN);
  mpz_init(biR);
  mpz_init(biS);
  mpz_init(biT);
  mpz_init(biU);
  mpz_init(biExp);
  mpz_init(biTmp);
}

/* ============================================================================================== */

void free_vars()
{
  int i = 0;
  for (i = 0 ; i < PWmax; i++)
  {
    mpz_clear(aiJS[i]);
    mpz_clear(aiJW[i]);
    mpz_clear(aiJX[i]);
    mpz_clear(aiJ0[i]);
    mpz_clear(aiJ1[i]);
    mpz_clear(aiJ2[i]);
    mpz_clear(aiJ00[i]);
    mpz_clear(aiJ01[i]);
  }
  free(aiJS);
  free(aiJW);
  free(aiJX);
  free(aiJ0);
  free(aiJ1);
  free(aiJ2);
  free(aiJ00);
  free(aiJ01);

  mpz_clear(TestNbr);
  mpz_clear(biN);
  mpz_clear(biR);
  mpz_clear(biS);
  mpz_clear(biT);
  mpz_clear(biU);
  mpz_clear(biExp);
  mpz_clear(biTmp);
}

/* ============================================================================================== */

// Compare Nbr1^2 vs. Nbr2
// -1 -> Nbr1^2 < Nbr2
//  0 -> Nbr1^2 == Nbr2
//  1 -> Nbr1^2 > Nbr2
int CompareSquare(mpz_t Nbr1, mpz_t Nbr2)
{
  mpz_t tmp;
  int cmp = 0;

  mpz_init(tmp);
  mpz_mul(tmp, Nbr1, Nbr1);

  cmp = mpz_cmp(tmp, Nbr2);
  mpz_clear(tmp);

  return cmp;
}

/* ============================================================================================== */

// Normalize coefficient of JS
void NormalizeJS(int PK, int PL, int PM, int P)
{
  int I, J;
  for (I = PL; I < PK; I++)
  {
    if (mpz_cmp_ui(aiJS[I], 0) != 0) /* (!BigNbrIsZero(aiJS[I])) */
    {
      /* biT = aiJS[I]; */
      mpz_set(biT, aiJS[I]);
      for (J = 1; J < P; J++)
      {
        /* SubtractBigNbrModN(aiJS[I - J * PM], biT, aiJS[I - J * PM], TestNbr, NumberLength); */
        mpz_sub(aiJS[I - J * PM], aiJS[I - J * PM], biT);
      }
      /* aiJS[I] = 0; */
      mpz_set_ui(aiJS[I], 0);
    }
  }
  for (I = 0; I < PK; I++)
    mpz_mod(aiJS[I], aiJS[I], TestNbr);
}

/* ============================================================================================== */

// Normalize coefficient of JW
void NormalizeJW(int PK, int PL, int PM, int P)
{
  int I, J;
  for (I = PL; I < PK; I++)
  {
    if (mpz_cmp_ui(aiJW[I], 0) != 0) /* (!BigNbrIsZero(aiJW[I])) */
    {
      /* biT = aiJW[I]; */
      mpz_set(biT, aiJW[I]);

      for (J = 1; J < P; J++)
      {
        /* SubtractBigNbrModN(aiJW[I - J * PM], biT, aiJW[I - J * PM], TestNbr, NumberLength); */
        mpz_sub(aiJW[I - J * PM], aiJW[I - J * PM], biT);
      }
      /* aiJW[I] = 0; */
      mpz_set_ui(aiJW[I], 0);
    }
  }
  for (I = 0; I < PK; I++)
    mpz_mod(aiJW[I], aiJW[I], TestNbr);
}

/* ============================================================================================== */

// Perform JS <- JS * JW

void JS_JW(int PK, int PL, int PM, int P)
{
  int I, J, K;
  for (I = 0; I < PL; I++)
  {
    for (J = 0; J < PL; J++)
    {
      K = (I + J) % PK;
      /* MontgomeryMult(aiJS[I], aiJW[J], biTmp); */
      /* AddBigNbrModN(aiJX[K], biTmp, aiJX[K], TestNbr, NumberLength); */
      mpz_mul(biTmp, aiJS[I], aiJW[J]);
      mpz_add(aiJX[K], aiJX[K], biTmp);
    }
  }
  for (I = 0; I < PK; I++)
  {
    /* aiJS[I] = aiJX[I]; */
    /* aiJX[I] = 0; */
    mpz_swap(aiJS[I], aiJX[I]);
    mpz_set_ui(aiJX[I], 0);
  }
  NormalizeJS(PK, PL, PM, P);
}

/* ============================================================================================== */

// Perform JS <- JS ^ 2

void JS_2(int PK, int PL, int PM, int P)
{
  int I, J, K;
  for (I = 0; I < PL; I++)
  {
    K = 2 * I % PK;
    /* MontgomeryMult(aiJS[I], aiJS[I], biTmp); */
    /* AddBigNbrModN(aiJX[K], biTmp, aiJX[K], TestNbr, NumberLength); */
    /* AddBigNbrModN(aiJS[I], aiJS[I], biT, TestNbr, NumberLength); */
    mpz_mul(biTmp, aiJS[I], aiJS[I]);
    mpz_add(aiJX[K], aiJX[K], biTmp);
    mpz_add(biT, aiJS[I], aiJS[I]);
    for (J = I + 1; J < PL; J++)
    {
      K = (I + J) % PK;
      /* MontgomeryMult(biT, aiJS[J], biTmp); */
      /* AddBigNbrModN(aiJX[K], biTmp, aiJX[K], TestNbr, NumberLength); */
      mpz_mul(biTmp, biT, aiJS[J]);
      mpz_add(aiJX[K], aiJX[K], biTmp);
    }
  }
  for (I = 0; I < PK; I++)
  {
    /* aiJS[I] = aiJX[I]; */
    /* aiJX[I] = 0; */
    mpz_swap(aiJS[I], aiJX[I]);
    mpz_set_ui(aiJX[I], 0);
  }
  NormalizeJS(PK, PL, PM, P);
}

/* ============================================================================================== */

// Perform JS <- JS ^ E

void JS_E(int PK, int PL, int PM, int P)
{
  int K;
  long Mask;

  if (mpz_cmp_ui(biExp, 1) == 0)
  {
    return;
  } // Return if E == 1


  for (K = 0; K < PL; K++)
  {
    mpz_set(aiJW[K], aiJS[K]);
  }


  Mask = mpz_sizeinbase(biExp, 2)-1;

  do
  {
    JS_2(PK, PL, PM, P);
    Mask--;
    if (mpz_tstbit(biExp, Mask))
    {
      JS_JW(PK, PL, PM, P);
    }
  }
  while (Mask > 0);
}

/* ============================================================================================== */

// if mode==0 then store J(p,q) ie x+f(x)
// if mode==1 then p=2, look for p==1, and stores Jstar(q) ie 2x+f(x)
// if mode==2 then p=2, look for p==4, and stores Jhash(q) ie 3x+f(x)
// This is based on ideas and code from Jason Moxham

void JacobiSum(int mode, int P, int PL, int Q)
{
  int I, a, myP;

  for (I = 0; I < PL; I++)
    mpz_set_ui(aiJ0[I], 0);

  myP = P; /* if (mode == 0) */
  if (mode == 1) myP = 1;
  if (mode == 2) myP = 4;

  for(a = 0; a < JPQSMAX; a++)
    if(jpqs[a].p == myP && jpqs[a].q == Q)
      break;

  for (I = 0; I < PL; I++)
    mpz_set_si(aiJ0[I], sls[jpqs[a].index+I]);
}
/* ============================================================================================== */

int mpz_aprcl(mpz_t N)
{
  /* Just return the status of the input */
  /* Do not print out any progress information */
  return mpz_aprtcle(N, APRTCLE_VERBOSE0);
}

/* ============================================================================================== */

/* Prime checking routine                       */
/* Return codes: 0 = N is composite.            */
/*               1 = N is a bpsw probable prime */
/*               2 = N is prime.                */
int mpz_aprtcle(mpz_t N, int verbose)
{
  int T, U;
  int i, j, H, I, J, K, P, Q, W, X;
  int IV, InvX, LEVELnow, NP, PK, PL, PM, SW, VK, TestedQs, TestingQs;
  int QQ, T1, T3, U1, U3, V1, V3;
  int break_this = 0;


  /* make sure the input is >= 2 and odd */
  if (mpz_cmp_ui(N, 2) < 0)
    return APRTCLE_COMPOSITE;

  if (mpz_divisible_ui_p(N, 2))
  {
    if (mpz_cmp_ui(N, 2) == 0)
      return APRTCLE_PRIME;
    else
      return APRTCLE_COMPOSITE;
  }

/* only three small exceptions for this implementation */
/* with this set of P and Q primes */
  if (mpz_cmp_ui(N, 3) == 0)
    return APRTCLE_PRIME;
  if (mpz_cmp_ui(N, 7) == 0)
    return APRTCLE_PRIME;
  if (mpz_cmp_ui(N, 11) == 0)
    return APRTCLE_PRIME;


  /* If the input number is larger than 4000 decimal digits
     we will just return whether it is a BPSW (probable) prime */
  NumberLength = mpz_sizeinbase(N, 10);
  if (NumberLength > 4000)
  {
    if (verbose >= APRTCLE_VERBOSE2)
      printf(" Info: Number too large, returning BPSW(N)\n");
    return mpz_bpsw_prp(N);
  }

  allocate_vars();

  mpz_set(TestNbr, N);
  mpz_set_si(biS, 0);

  j = PK = PL = PM = 0;
  for (J = 0; J < PWmax; J++)
  {
    /* aiJX[J] = 0; */
    mpz_set_ui(aiJX[J], 0);
  }
  break_this = 0;
/* GetPrimes2Test : */
  for (i = 0; i < LEVELmax; i++)
  {
    /* biS[0] = 2; */
    mpz_set_ui(biS, 2);
    
    for (j = 0; j < aiNQ[i]; j++)
    {
      Q = aiQ[j];
      if (aiT[i]%(Q-1) != 0) continue;
/*
      U = aiT[i] * Q;
      do
      {
        U /= Q;
        // MultBigNbrByLong(biS, Q, biS, NumberLength);
        mpz_mul_ui(biS, biS, Q);
      }
      while (U % Q == 0);
*/
      mpz_set_ui(biU, aiT[i]);
      mpz_mul_ui(biU, biU, Q);
      do
      {
        mpz_divexact_ui(biU, biU, Q);
        // MultBigNbrByLong(biS, Q, biS, NumberLength);
        mpz_mul_ui(biS, biS, Q);
      }
      while (mpz_divisible_ui_p(biU, Q));

      // Exit loop if S^2 > N.
      if (CompareSquare(biS, TestNbr) > 0)
      {
        /* break GetPrimes2Test; */
        break_this = 1;
        break;
      }
    } /* End for j */
    if (break_this) break;
  } /* End for i */
  if (i == LEVELmax)
  { /* too big */
    free_vars();
    return mpz_bpsw_prp(N);
  }
  LEVELnow = i;
  TestingQs = j;
  T = aiT[LEVELnow];
  NP = aiNP[LEVELnow];

MainStart:
  for (;;)
  {
    for (i = 0; i < NP; i++)
    {
      P = aiP[i];
      if (T%P != 0) continue;

      SW = TestedQs = 0;
      /* Q = W = (int) BigNbrModLong(TestNbr, P * P); */
      Q = W = mpz_fdiv_ui(TestNbr, P * P);
      for (J = P - 2; J > 0; J--)
      {
        W = (W * Q) % (P * P);
      }
      if (P > 2 && W != 1)
      {
        SW = 1;
      }
      for (;;)
      {
        for (j = TestedQs; j <= TestingQs; j++)
        {
          Q = aiQ[j] - 1;
          K = 0;
          while (Q % P == 0)
          {
            K++;
            Q /= P;
          }
          Q = aiQ[j];
          if (K == 0)
          {
            continue;
          }

          if (verbose >= APRTCLE_VERBOSE1)
          {
            printf("P = %2d, Q = %12d  (%3.2f%%)\r", P, Q, (i * (TestingQs + 1) + j) * 100.0 / (NP * (TestingQs + 1)));
            fflush(stdout);
          }

          PM = 1;
          for (I = 1; I < K; I++)
          {
            PM = PM * P;
          }
          PL = (P - 1) * PM;
          PK = P * PM;
          for (I = 0; I < PK; I++)
          {
            /* aiJ0[I] = aiJ1[I] = 0; */
            mpz_set_ui(aiJ0[I], 0);
            mpz_set_ui(aiJ1[I], 0);
          }
          if (P > 2)
          {
            JacobiSum(0, P, PL, Q);
          }
          else
          {
            if (K != 1)
            {
              JacobiSum(0, P, PL, Q);
              for (I = 0; I < PK; I++)
              {
                /* aiJW[I] = 0; */
                mpz_set_ui(aiJW[I], 0);
              }
              if (K != 2)
              {
                for (I = 0; I < PM; I++)
                {
                  /* aiJW[I] = aiJ0[I]; */
                  mpz_set(aiJW[I], aiJ0[I]);
                }
                JacobiSum(1, P, PL, Q);
                for (I = 0; I < PM; I++)
                {
                  /* aiJS[I] = aiJ0[I]; */
                  mpz_set(aiJS[I], aiJ0[I]);
                }
                JS_JW(PK, PL, PM, P);
                for (I = 0; I < PM; I++)
                {
                  /* aiJ1[I] = aiJS[I]; */
                  mpz_set(aiJ1[I], aiJS[I]);
                }
                JacobiSum(2, P, PL, Q);
                for (I = 0; I < PK; I++)
                {
                  /* aiJW[I] = 0; */
                  mpz_set_ui(aiJW[I], 0);
                }
                for (I = 0; I < PM; I++)
                {
                  /* aiJS[I] = aiJ0[I]; */
                  mpz_set(aiJS[I], aiJ0[I]);
                }
                JS_2(PK, PL, PM, P);
                for (I = 0; I < PM; I++)
                {
                  /* aiJ2[I] = aiJS[I]; */
                  mpz_set(aiJ2[I], aiJS[I]);
                }
              }
            }
          }
          /* aiJ00[0] = aiJ01[0] = 1; */
          mpz_set_ui(aiJ00[0], 1);
          mpz_set_ui(aiJ01[0], 1);
          for (I = 1; I < PK; I++)
          {
            /* aiJ00[I] = aiJ01[I] = 0; */
            mpz_set_ui(aiJ00[I], 0);
            mpz_set_ui(aiJ01[I], 0);
          }
          /* VK = (int) BigNbrModLong(TestNbr, PK); */
          VK = mpz_fdiv_ui(TestNbr, PK);
          for (I = 1; I < PK; I++)
          {
            if (I % P != 0)
            {
              U1 = 1;
              U3 = I;
              V1 = 0;
              V3 = PK;
              while (V3 != 0)
              {
                QQ = U3 / V3;
                T1 = U1 - V1 * QQ;
                T3 = U3 - V3 * QQ;
                U1 = V1;
                U3 = V3;
                V1 = T1;
                V3 = T3;
              }
              aiInv[I] = (U1 + PK) % PK;
            }
            else
            {
              aiInv[I] = 0;
            }
          }
          if (P != 2)
          {
            for (IV = 0; IV <= 1; IV++)
            {
              for (X = 1; X < PK; X++)
              {
                for (I = 0; I < PK; I++)
                {
                  /* aiJS[I] = aiJ0[I]; */
                  mpz_set(aiJS[I], aiJ0[I]);
                }
                if (X % P == 0)
                {
                  continue;
                }
                if (IV == 0)
                {
                  /* LongToBigNbr(X, biExp, NumberLength); */
                  mpz_set_ui(biExp, X);
                }
                else
                {
                  /* LongToBigNbr(VK * X / PK, biExp, NumberLength); */
                  mpz_set_ui(biExp, (VK * X) / PK);
                  if ((VK * X) / PK == 0)
                  {
                    continue;
                  }
                }
                JS_E(PK, PL, PM, P);
                for (I = 0; I < PK; I++)
                {
                  /* aiJW[I] = 0; */
                  mpz_set_ui(aiJW[I], 0);
                }
                InvX = aiInv[X];
                for (I = 0; I < PK; I++)
                {
                  J = (I * InvX) % PK;
                  /* AddBigNbrModN(aiJW[J], aiJS[I], aiJW[J], TestNbr, NumberLength); */
                  mpz_add(aiJW[J], aiJW[J], aiJS[I]);
                }
                NormalizeJW(PK, PL, PM, P);
                if (IV == 0)
                {
                  for (I = 0; I < PK; I++)
                  {
                    /* aiJS[I] = aiJ00[I]; */
                    mpz_set(aiJS[I], aiJ00[I]);
                  }
                }
                else
                {
                  for (I = 0; I < PK; I++)
                  {
                    /* aiJS[I] = aiJ01[I]; */
                    mpz_set(aiJS[I], aiJ01[I]);
                  }
                }
                JS_JW(PK, PL, PM, P);
                if (IV == 0)
                {
                  for (I = 0; I < PK; I++)
                  {
                    /* aiJ00[I] = aiJS[I]; */
                    mpz_set(aiJ00[I], aiJS[I]);
                  }
                }
                else
                {
                  for (I = 0; I < PK; I++)
                  {
                    /* aiJ01[I] = aiJS[I]; */
                    mpz_set(aiJ01[I], aiJS[I]);
                  }
                }
              } /* end for X */
            } /* end for IV */
          }
          else
          {
            if (K == 1)
            {
              /* MultBigNbrByLongModN(1, Q, aiJ00[0], TestNbr, NumberLength); */
              mpz_set_ui(aiJ00[0], Q);
              /* aiJ01[0] = 1; */
              mpz_set_ui(aiJ01[0], 1);
            }
            else
            {
              if (K == 2)
              {
                if (VK == 1)
                {
                  /* aiJ01[0] = 1; */
                  mpz_set_ui(aiJ01[0], 1);
                }
                /* aiJS[0] = aiJ0[0]; */
                /* aiJS[1] = aiJ0[1]; */
                mpz_set(aiJS[0], aiJ0[0]);
                mpz_set(aiJS[1], aiJ0[1]);
                JS_2(PK, PL, PM, P);
                if (VK == 3)
                {
                  /* aiJ01[0] = aiJS[0]; */
                  /* aiJ01[1] = aiJS[1]; */
                  mpz_set(aiJ01[0], aiJS[0]);
                  mpz_set(aiJ01[1], aiJS[1]);
                }
                /* MultBigNbrByLongModN(aiJS[0], Q, aiJ00[0], TestNbr, NumberLength); */
                mpz_mul_ui(aiJ00[0], aiJS[0], Q);
                /* MultBigNbrByLongModN(aiJS[1], Q, aiJ00[1], TestNbr, NumberLength); */
                mpz_mul_ui(aiJ00[1], aiJS[1], Q);
              }
              else
              {
                for (IV = 0; IV <= 1; IV++)
                {
                  for (X = 1; X < PK; X += 2)
                  {
                    for (I = 0; I <= PM; I++)
                    {
                      /* aiJS[I] = aiJ1[I]; */
                      mpz_set(aiJS[I], aiJ1[I]);
                    }
                    if (X % 8 == 5 || X % 8 == 7)
                    {
                      continue;
                    }
                    if (IV == 0)
                    {
                      /* LongToBigNbr(X, biExp, NumberLength); */
                      mpz_set_ui(biExp, X);
                    }
                    else
                    {
                      /* LongToBigNbr(VK * X / PK, biExp, NumberLength); */
                      mpz_set_ui(biExp, VK * X / PK);
                      if (VK * X / PK == 0)
                      {
                        continue;
                      }
                    }
                    JS_E(PK, PL, PM, P);
                    for (I = 0; I < PK; I++)
                    {
                      /* aiJW[I] = 0; */
                      mpz_set_ui(aiJW[I], 0);
                    }
                    InvX = aiInv[X];
                    for (I = 0; I < PK; I++)
                    {
                      J = I * InvX % PK;
                      /* AddBigNbrModN(aiJW[J], aiJS[I], aiJW[J], TestNbr, NumberLength); */
                      mpz_add(aiJW[J], aiJW[J], aiJS[I]);
                    }
                    NormalizeJW(PK, PL, PM, P);
                    if (IV == 0)
                    {
                      for (I = 0; I < PK; I++)
                      {
                        /* aiJS[I] = aiJ00[I]; */
                        mpz_set(aiJS[I], aiJ00[I]);
                      }
                    }
                    else
                    {
                      for (I = 0; I < PK; I++)
                      {
                        /* aiJS[I] = aiJ01[I]; */
                        mpz_set(aiJS[I], aiJ01[I]);
                      }
                    }
                    NormalizeJS(PK, PL, PM, P);
                    JS_JW(PK, PL, PM, P);
                    if (IV == 0)
                    {
                      for (I = 0; I < PK; I++)
                      {
                        /* aiJ00[I] = aiJS[I]; */
                        mpz_set(aiJ00[I], aiJS[I]);
                      }
                    }
                    else
                    {
                      for (I = 0; I < PK; I++)
                      {
                        /* aiJ01[I] = aiJS[I]; */
                        mpz_set(aiJ01[I], aiJS[I]);
                      }
                    }
                  } /* end for X */
                  if (IV == 0 || VK % 8 == 1 || VK % 8 == 3)
                  {
                    continue;
                  }
                  for (I = 0; I < PM; I++)
                  {
                    /* aiJW[I] = aiJ2[I]; */
                    /* aiJS[I] = aiJ01[I]; */
                    mpz_set(aiJW[I], aiJ2[I]);
                    mpz_set(aiJS[I], aiJ01[I]);
                  }
                  for (; I < PK; I++)
                  {
                    /* aiJW[I] = aiJS[I] = 0; */
                    mpz_set_ui(aiJW[I], 0);
                    mpz_set_ui(aiJS[I], 0);
                  }
                  JS_JW(PK, PL, PM, P);
                  for (I = 0; I < PM; I++)
                  {
                    /* aiJ01[I] = aiJS[I]; */
                    mpz_set(aiJ01[I], aiJS[I]);
                  }
                } /* end for IV */
              }
            }
          }
          for (I = 0; I < PL; I++)
          {
            /* aiJS[I] = aiJ00[I]; */
            mpz_set(aiJS[I], aiJ00[I]);
          }
          for (; I < PK; I++)
          {
            /* aiJS[I] = 0; */
            mpz_set_ui(aiJS[I], 0);
          }
          /* DivBigNbrByLong(TestNbr, PK, biExp, NumberLength); */
          mpz_fdiv_q_ui(biExp, TestNbr, PK);
          JS_E(PK, PL, PM, P);
          for (I = 0; I < PK; I++)
          {
            /* aiJW[I] = 0; */
            mpz_set_ui(aiJW[I], 0);
          }
          for (I = 0; I < PL; I++)
          {
            for (J = 0; J < PL; J++)
            {
              /* MontgomeryMult(aiJS[I], aiJ01[J], biTmp); */
              /* AddBigNbrModN(biTmp, aiJW[(I + J) % PK], aiJW[(I + J) % PK], TestNbr, NumberLength); */
              mpz_mul(biTmp, aiJS[I], aiJ01[J]);
              mpz_add(aiJW[(I + J) % PK], biTmp, aiJW[(I + J) % PK]);
            }
          }
          NormalizeJW(PK, PL, PM, P);
/* MatchingRoot : */
          do
          {
            H = -1;
            W = 0;
            for (I = 0; I < PL; I++)
            {
              if (mpz_cmp_ui(aiJW[I], 0) != 0)/* (!BigNbrIsZero(aiJW[I])) */
              {
                /* if (H == -1 && BigNbrAreEqual(aiJW[I], 1)) */
                if (H == -1 && (mpz_cmp_ui(aiJW[I], 1) == 0))
                {
                  H = I;
                }
                else
                {
                  H = -2;
                  /* AddBigNbrModN(aiJW[I], MontgomeryMultR1, biTmp, TestNbr, NumberLength); */
                  mpz_add_ui(biTmp, aiJW[I], 1);
                  mpz_mod(biTmp, biTmp, TestNbr);
                  if (mpz_cmp_ui(biTmp, 0) == 0) /* (BigNbrIsZero(biTmp)) */
                  {
                    W++;
                  }
                }
              }
            }
            if (H >= 0)
            {
              /* break MatchingRoot; */
              break;
            }
            if (W != P - 1)
            {
              /* Not prime */
              if (verbose >= APRTCLE_VERBOSE2)
                {printf("Failed: W != P - 1 : H=%d : W=%d : P-1=%d\n", H, W, P-1); fflush(stdout);}
              free_vars();
              return APRTCLE_COMPOSITE;
            }
            for (I = 0; I < PM; I++)
            {
              /* AddBigNbrModN(aiJW[I], 1, biTmp, TestNbr, NumberLength); */
              mpz_add_ui(biTmp, aiJW[I], 1);
              mpz_mod(biTmp, biTmp, TestNbr);
              if (mpz_cmp_ui(biTmp, 0) == 0) /* (BigNbrIsZero(biTmp)) */
              {
                break;
              }
            }
            if (I == PM)
            {
              /* Not prime */
              if (verbose >= APRTCLE_VERBOSE2)
                {printf("Failed: I == PM : I=%d : PM=%d\n", I, PM); fflush(stdout);}
              free_vars();
              return APRTCLE_COMPOSITE;
            }
            for (J = 1; J <= P - 2; J++)
            {
              /* AddBigNbrModN(aiJW[I + J * PM], 1, biTmp, TestNbr, NumberLength); */
              mpz_add_ui(biTmp, aiJW[I + J * PM], 1);
              mpz_mod(biTmp, biTmp, TestNbr);
              if (mpz_cmp_ui(biTmp, 0) != 0)/* (!BigNbrIsZero(biTmp)) */
              {
                /* Not prime */
                if (verbose >= APRTCLE_VERBOSE2)
                {
                  printf("Failed: biTmp != 0 (1)\n");
                  gmp_printf("biTmp=%Zd\n", biTmp); fflush(stdout);
                }
                free_vars();
                return APRTCLE_COMPOSITE;
              }
            }
            H = I + PL;
          }
          while (0);

          if (SW == 1 || H % P == 0)
          {
            continue;
          }
          if (P != 2)
          {
            SW = 1;
            continue;
          }
          if (K == 1)
          {
            if ((mpz_get_ui(TestNbr) & 3) == 1)
            {
              SW = 1;
            }
            continue;
          }

          // if (Q^((N-1)/2) mod N != N-1), N is not prime.

          /* MultBigNbrByLongModN(1, Q, biTmp, TestNbr, NumberLength); */
          mpz_set_ui(biTmp, Q);
          mpz_mod(biTmp, biTmp, TestNbr);

          mpz_sub_ui(biT, TestNbr, 1); /* biT = n-1 */
          mpz_divexact_ui(biT, biT, 2); /* biT = (n-1)/2 */
          mpz_powm(biR, biTmp, biT, TestNbr); /* biR = Q^((n-1)/2) mod n */
          mpz_add_ui(biTmp, biR, 1);
          mpz_mod(biTmp, biTmp, TestNbr);

          if (mpz_cmp_ui(biTmp, 0) != 0)/* (!BigNbrIsZero(biTmp)) */
          {
            /* Not prime */
            if (verbose >= APRTCLE_VERBOSE2)
            {
              printf("Failed: biTmp != 0 (2)\n");
              gmp_printf("biTmp = %Zd\n", biTmp); fflush(stdout);
            }
            free_vars();
            return APRTCLE_COMPOSITE;
          }
          SW = 1;
        } /* end for j */
        if (SW == 0)
        {
          TestedQs = TestingQs + 1;
          if (TestingQs < aiNQ[LEVELnow] - 1)
          {
            TestingQs++;
            Q = aiQ[TestingQs];
/*
            U = T * Q;
            do
            {
              // MultBigNbrByLong(biS, Q, biS, NumberLength);
              mpz_mul_ui(biS, biS, Q);
              U /= Q;
            }
            while (U % Q == 0);
*/
            mpz_set_ui(biU, T);
            mpz_mul_ui(biU, biU, Q);
            do
            {
              // MultBigNbrByLong(biS, Q, biS, NumberLength);
              mpz_mul_ui(biS, biS, Q);
              mpz_divexact_ui(biU, biU, Q);
            }
            while (mpz_divisible_ui_p(biU, Q));

            continue; /* Retry */
          }
          LEVELnow++;
          if (LEVELnow == LEVELmax)
          {
            if (verbose >= APRTCLE_VERBOSE2)
              {printf("APR-CL cannot tell: lvlnow == lvlmax, returning BPSW(N)\n"); fflush(stdout);}
            free_vars();
            return mpz_bpsw_prp(N); /* Cannot tell */
          }
          T = aiT[LEVELnow];
          NP = aiNP[LEVELnow];
          /* biS = 2; */
          mpz_set_ui(biS, 2);
          for (J = 0; J <= aiNQ[LEVELnow]; J++)
          {
            Q = aiQ[J];
            if (T%(Q-1) != 0) continue;
/*
            U = T * Q;
            do
            {
              // MultBigNbrByLong(biS, Q, biS, NumberLength);
              mpz_mul_ui(biS, biS, Q);
              U /= Q;
            }
            while (U % Q == 0);
*/
            mpz_set_ui(biU, T);
            mpz_mul_ui(biU, biU, Q);
            do
            {
              // MultBigNbrByLong(biS, Q, biS, NumberLength);
              mpz_mul_ui(biS, biS, Q);
              mpz_divexact_ui(biU, biU, Q);
            }
            while (mpz_divisible_ui_p(biU, Q));

            if (CompareSquare(biS, TestNbr) > 0)
            {
              TestingQs = J;
              /* continue MainStart; */ /* Retry from the beginning */
              goto MainStart;
            }
          } /* end for J */
          if (verbose >= APRTCLE_VERBOSE2)
            {printf("Failed: APR-CL error, returning BPSW(N)\n"); fflush(stdout);}
          free_vars();
          return mpz_bpsw_prp(N); /* Program error */
        } /* end if */
        break;
      } /* end for (;;) */
    } /* end for i */

    // Final Test
   
    /* biR = 1 */
    mpz_set_ui(biR, 1);
    /* biN <- TestNbr mod biS */ /* Compute N mod S */
    mpz_fdiv_r(biN, TestNbr, biS);
   
    for (U = 1; U <= T; U++)
    {
      /* biR <- (biN * biR) mod biS */
      mpz_mul(biR, biN, biR);
      mpz_mod(biR, biR, biS);
      if (mpz_cmp_ui(biR, 1) == 0) /* biR == 1 */
      {
        /* Number is prime */
        free_vars();
        return APRTCLE_PRIME;
      }
      if (mpz_divisible_p(TestNbr, biR) && mpz_cmp(biR, TestNbr) < 0) /* biR < N and biR | TestNbr */
      {
        /* Number is composite */
        if (verbose >= APRTCLE_VERBOSE2)
          {gmp_printf(" *** Found factor: %Zd\n", biR); fflush(stdout);}
        free_vars();
        return APRTCLE_COMPOSITE;
      }
    } /* End for U */
    /* This should never be reached. */
    if (verbose >= APRTCLE_VERBOSE2)
      {printf("Failed: APR-CL error with final test, returning BPSW(N)\n"); fflush(stdout);}
    free_vars();
    return mpz_bpsw_prp(N); /* Program error */
  }
}

/* ============================================================================================== */

// Minimal public MPFR ABI declarations for systems with the runtime but no headers.
// Prefer the installed development header. Fallback supports LP64 Unix only.
#pragma once
#if __has_include(<mpfr.h>)
#include <mpfr.h>
#else
#include <climits>
static_assert(sizeof(long)==8 && sizeof(unsigned long)==8, "Install MPFR headers on non-LP64 platforms");
extern "C" {
typedef long mpfr_prec_t;
typedef long mpfr_exp_t;
typedef struct { mpfr_prec_t _mpfr_prec; int _mpfr_sign; mpfr_exp_t _mpfr_exp; unsigned long *_mpfr_d; } __mpfr_struct;
typedef __mpfr_struct mpfr_t[1];
typedef __mpfr_struct *mpfr_ptr;
typedef const __mpfr_struct *mpfr_srcptr;
typedef enum { MPFR_RNDN=0, MPFR_RNDZ, MPFR_RNDU, MPFR_RNDD, MPFR_RNDA, MPFR_RNDF } mpfr_rnd_t;
void mpfr_init2(mpfr_ptr,mpfr_prec_t);
void mpfr_clear(mpfr_ptr);
int mpfr_set(mpfr_ptr,mpfr_srcptr,mpfr_rnd_t);
int mpfr_set_si(mpfr_ptr,long,mpfr_rnd_t);
int mpfr_set_str(mpfr_ptr,const char*,int,mpfr_rnd_t);
int mpfr_add(mpfr_ptr,mpfr_srcptr,mpfr_srcptr,mpfr_rnd_t);
int mpfr_sub(mpfr_ptr,mpfr_srcptr,mpfr_srcptr,mpfr_rnd_t);
int mpfr_mul(mpfr_ptr,mpfr_srcptr,mpfr_srcptr,mpfr_rnd_t);
int mpfr_div(mpfr_ptr,mpfr_srcptr,mpfr_srcptr,mpfr_rnd_t);
int mpfr_neg(mpfr_ptr,mpfr_srcptr,mpfr_rnd_t);
int mpfr_abs(mpfr_ptr,mpfr_srcptr,mpfr_rnd_t);
int mpfr_sqrt(mpfr_ptr,mpfr_srcptr,mpfr_rnd_t);
int mpfr_cmp(mpfr_srcptr,mpfr_srcptr);
double mpfr_get_d(mpfr_srcptr,mpfr_rnd_t);
int mpfr_snprintf(char*,size_t,const char*,...);
const char *mpfr_get_version(void);
}
#endif

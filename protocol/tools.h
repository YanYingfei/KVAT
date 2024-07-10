#ifndef TOOLS_H
#define TOOLS_H
#include <stdint.h>
#include "../poly/polymat.h"
#include "../poly/poly.h"
#include "../poly/fips202.h"
#include "params.h"
void randseed(unsigned char* seed , int seedlen);
void matrix_mul(int64_t *res , int64_t *matrix , int64_t *vector , int k , int l);
void vec_add(int64_t *a , int64_t *b , int len);
void vec_mul_si(int64_t *a , int64_t x , int len);
int64_t vec_mul(int64_t *a , int64_t *b , int len);
int64_t vec_norm(int64_t *vec , int veclen);
int Rej0(int64_t *z , int64_t *v , double s , double GAMMA);
int Rej1(polyvecq *z ,polyvecq *v ,double s , double GAMMA);
void sample_gamma(keccak_state *state , polymatq *res);
void vec2char(unsigned char* res , int64_t *vec , int len);
void char2vec(int64_t *res , unsigned char *vecchar, int len);
void polyvec2vec(int64_t *res  , polyvecq *y);
void geta(Polyq *ai , int norm);
int Cjudge(Polyq* a);
void Csamplepoly(keccak_state*state , Polyq *res);
#endif
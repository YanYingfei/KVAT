#ifndef SAMPLE_H
#define SAMPLE_H
#include "fips202.h"
#include "params.h"
void uniformq(keccak_state *state , int64_t *res , int len = 128);
void uniformp(keccak_state* state, int16_t *res);

template <typename T>
void cbd(keccak_state*state , T *res , int eta ,int len=128);
template <typename T>
void Gaussian(T *res , double sigma);
#endif
#ifndef PROTOCOL_batch_H
#define PROTOCOL_batch_H
#include "../poly/polymat.h"
#include "params.h"

int Prove_batch(unsigned char* pi , unsigned char* crsseed , int b[],polyvecq *s0, polyvecq *s1,polyvecq *e0,  polyvecq *e1, polyvecq *e3, \
polyvecq* h, polymatq*A, polyvecq*t0, polyvecq *t1, polymatq *c);
int Verify_batch(unsigned char *crsseed, unsigned char *pi , polymatq *A, polymatq *c , polyvecq *t0 , polyvecq* t1 , polyvecq *h);
#endif
#ifndef PROTOCOL2_H
#define PROTOCOL2_H
#include "../poly/polymat.h"
#include "params.h"

int Prove2(unsigned char* pi , unsigned char* crsseed ,polyvecq *sb , polyvecq *s2 , polyvecq *e1 , Polyq *e3, Polyq* h,polymatq*A,polyvecq*t, polyvecq *b,\
    polyvecq *c);
int Verify2(unsigned char *crsseed, unsigned char *pi , polymatq *A, polyvecq *b, polyvecq *c , polyvecq *t , Polyq *h);
#endif
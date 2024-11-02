#ifndef PROTOCOL2_H
#define PROTOCOL2_H
#include "../poly/polymat.h"
#include "params.h"

int Prove2(unsigned char* pi , unsigned char* crsseed , int b ,polyvecq *sb , polyvecq *e1 , Polyq *e3, Polyq* h,polymatq*A,polyvecq*t0, polyvecq *t1,\
    polyvecq *c);
int Verify2(unsigned char *crsseed, unsigned char *pi , polymatq *A, polyvecq *c , polyvecq *t0 , polyvecq* t1 , Polyq *h);
#endif
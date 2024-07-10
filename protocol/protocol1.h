#ifndef PROTOCOL1_H
#define PROTOCOL1_H
#include "../poly/polymat.h"
#include "params.h"

int Prove1(unsigned char* pi , unsigned char* crsseed ,unsigned char*rseed,polyvecq *x , polyvecq* e2, unsigned char* H1input, unsigned char* uv ,\
    unsigned char*pke , polymatq *A , polyvecq *c);
void setup(unsigned char* crsseed);
int Verify1(unsigned char *crsseed, unsigned char *pi , polymatq *A, polyvecq *c);
#endif
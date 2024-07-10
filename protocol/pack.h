#ifndef PACK_H
#define PACK_H

#include "../poly/poly.h"
#include "../poly/polymat.h"
#include "../poly/params.h"
#include "params.h"
#include <cstring>

void packpk(unsigned char* pk ,unsigned char* seedA, unsigned char *seedb , polyvecq*t);
void unpackpk(unsigned char* pk , unsigned char* seedA , unsigned char* seedb, polyvecq*t);
void packsk(unsigned char*sk ,polyvecq*s0,polyvecq*s1,polyvecq*s2,polyvecq*e1 , unsigned char* k);
void unpacksk(unsigned char*sk ,polyvecq*s0,polyvecq*s1,polyvecq*s2,polyvecq*e1 , unsigned char* k);
void Compress(unsigned char*res, Polyq* h );
void Decompress(Polyq* h , unsigned char*hchar);
#endif
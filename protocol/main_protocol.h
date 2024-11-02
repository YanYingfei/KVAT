#ifndef MAIN_PROTOCOL_H
#define MAIN_PROTOCOL_H

#include "../poly/poly.h"
#include "../poly/polymat.h"
#include "../poly/params.h"
#include "protocol1.h"
#include "protocol2.h"
#include "cpaenc.h"
#include "tools.h"
#include "../poly/fips202.h"


void keygen(unsigned char* pk , unsigned char*sk);
void AT_CQ(unsigned char* pi1 ,unsigned char*crsseed, unsigned char* query, unsigned char* H1input ,unsigned char* xchar , unsigned char* pk , double *pitime = NULL);
void AT_IT(unsigned char* pi2 , unsigned char* crsseed2 , unsigned char* resp ,unsigned char* pk, unsigned char *sk ,\
     unsigned char* query , unsigned char* crsseed , unsigned char* pi1 , double *pi1time = NULL , double *pi2time = NULL);
void AT_CF(unsigned char* msg , unsigned char* pk , unsigned char*resp ,unsigned char*H1input, unsigned char* xchar ,\
    unsigned char* query , unsigned char* pi2 , unsigned char* crsseed2 , double *pitime = NULL);
int AT_Rb(unsigned char* pk , unsigned char* sk , unsigned char* msg);

#endif
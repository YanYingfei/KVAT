#include "../poly/fips202.h"
#include "../poly/poly.h"
#include "../poly/polymat.h"
#include "params.h"
#include <ctime>
#include <stdlib.h>
#include "cpaenc.h"
#include <cmath>
#include "tools.h"
void keygen(unsigned char* pke){
    randseed(pke , 64);
}
void CompressPoly(unsigned char*res , Polyp *x , int d){
    int dmod;
    int32_t temp[128];
    if(d == 10){
        dmod = 1024;
    }
    if(d == 4){
        dmod = 16;
    }
    for(int j = 0 ; j < 128 ; j++){
        temp[j] = (int)round((double)((int32_t)x->polyarray[j]*dmod)/3329);
        if(temp[j] == dmod){
            temp[j] = 0;
        }
    }
    if(d == 4){
        for(int i = 0 ; i < 64 ; i++){
            res[i] = (temp[i*2]<<4)|(temp[i*2+1]);
        }
    }
    else{
        for(int i = 0 ; i < 32;i++){
            res[i*5] = (temp[4*i]>>2);
            res[i*5+1] = ((temp[4*i]&3)<<6)|(temp[4*i+1]>>4);
            res[i*5+2] = ((temp[4*i+1]&0xf)<<4)|(temp[4*i+2]>>6);
            res[i*5+3] = ((temp[4*i+2]&0x3f)<<2)|(temp[4*i+3]>>8);
            res[i*5+4] = temp[4*i+3]&0xff; 
        }
    }
}
void encodem(Polyp* mpoly , unsigned char*m){
    for(int i = 0 ; i < 16;i++){
        for(int j = 7 ; j >=0;j--){
            if((m[i]>>j)&1 == 1){
                mpoly->polyarray[i*8+j] += 3329/2;
            }
        }
    }
}
void Compressvec(unsigned char*res , polyvecp *u){
    int32_t temp[128];
    
    for(int i = 0 ; i < 4; i++){
        CompressPoly(res+UPOLYLEN*i , u->polyarray[i] ,10);
    }

}
void enc(unsigned char *cipher,unsigned char* uv , unsigned char* rseed, unsigned char *pke, unsigned char*m , int mlen){
    int blocklen = mlen / 16;
    polymatp *A = new polymatp(pke,32,0,4,4,1);
    keccak_state state1 , state2;
    shake256_init(&state1);
    shake256_absorb(&state1 , pke+32 , 32);
    polyvecp *t;

    unsigned char seed[32];
    randseed(seed , 32);
    memcpy(rseed , seed , 32);
    shake256_init(&state2);
    shake256_absorb(&state2 , seed , 32);
    polyvecp *r = new polyvecp(&state2 , 3 , 4 , 0);
    polyvecp *e1 = new polyvecp(&state2 , 2 , 4 , 0);
    polyvecp *u = new polyvecp(4);
    A->right_mul(u , r);
    u->to_poly();
    u->add(u , e1);
    Compressvec(cipher , u);
    u->to_char(uv);
    Polyp *e2 , *temp = new Polyp(1);
    for(int i = 0;i < blocklen ;i++){
        t = new polyvecp(&state1 , 0 , 4 , 1);
        t->mul(temp ,r);
        e2 = new Polyp(&state2 , 2 , 0);
        temp->to_poly();
        temp->add(temp , e2);
        encodem(temp , m+i*16);
        CompressPoly(cipher+UVECLEN+i*VPOLYLEN , temp , 4);
        temp->to_char(uv+POLYPLEN*4+i*POLYPLEN , 0);
        delete t;
    }

}
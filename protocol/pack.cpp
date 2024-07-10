
#include "pack.h"
#include <cmath>
#include <iostream>

void packpk(unsigned char* pk ,unsigned char* seedA, unsigned char *seedb , polyvecq*t){
    memcpy(pk , seedA , 32);
    memcpy(pk+32 , seedb , 32);
    t->to_char(pk+64);
}

void unpackpk(unsigned char* pk , unsigned char* seedA , unsigned char* seedb, polyvecq*t){
    memcpy(seedA , pk , 32);
    memcpy(seedb , pk+32, 32);
    t->from_char(pk+64 , 1);
}

void packsk(unsigned char*sk ,polyvecq*s0,polyvecq*s1,polyvecq*s2,polyvecq*e1 ,unsigned char* k){
    s0->to_char(sk , ETAR);
    s1->to_char(sk+VECETALEN , ETAR);
    s2->to_char(sk+2*VECETALEN , ETAR);
    e1->to_char(sk+3*VECETALEN , ETAR);
    memcpy(sk+4*VECETALEN ,k , 32);
}

void unpacksk(unsigned char*sk ,polyvecq*s0,polyvecq*s1,polyvecq*s2,polyvecq*e1 , unsigned char* k){
    s0->from_char(sk , 0 , ETAR);
    s1->from_char(sk+VECETALEN ,0 ,  ETAR);
    s2->from_char(sk+2*VECETALEN ,0 ,  ETAR);
    e1->from_char(sk+3*VECETALEN ,0 ,  ETAR);
    memcpy(k , sk+4*VECETALEN , 32);
}

void Compress(unsigned char *res, Polyq *h){
    int temp[128];
    h->to_poly();
    for(int i = 0 ; i < 128 ;i++){
        temp[i] = (int)round((double)(h->mod(h->polyarray[i])*8)/Q);
        if(temp[i] == 8){
            temp[i] = 0;
        } 
    }
    for(int i = 0 ; i < 16;i++){
        res[i*3] = (temp[i*8] << 5) | (temp[i*8+1]<<2) | (temp[i*8+2]>>1);
        res[i*3+1] = ((temp[i*8+2]&1) << 7) | (temp[i*8+3] << 4) | (temp[i*8+4]<<1) | (temp[i*8+5] >>2);
        res[i*3+2] = ((temp[i*8+5]&3) << 6) | (temp[i*8+6] << 3) | temp[i*8+7];
    }
}

void Decompress(Polyq *h, unsigned char *hchar){
    int temp[8];
    for(int i = 0 ; i < 16;i++){
        temp[0] = (hchar[i*3] >> 5)&7;
        temp[1] = (hchar[i*3] >> 2)&7;
        temp[2] = ((hchar[i*3]&3)<<1) | ((hchar[i*3+1]>>7));
        temp[3] = (hchar[i*3+1] >> 4)&7;
        temp[4] = (hchar[i*3+1] >> 1)&7;
        temp[5] = ((hchar[i*3+1]&1)<<2) | (hchar[i*3+2]>>6);
        temp[6] = (hchar[i*3+2] >> 3)&7;
        temp[7] = (hchar[i*3+2])&7;
        for(int j = 0 ; j < 8 ; j++){
            h->polyarray[i*8+j] = h->mod(Q*(temp[j])/8);
            //h->polyarray[i*8+j] =temp[j]-4;
        }
    }
    h->nttflag = 0;
}

#include "../poly/poly.h"
#include "../poly/polymat.h"
#include "../poly/params.h"
#include "pack.h"
#include "params.h"
#include "../poly/fips202.h"
#include <stdlib.h>
#include <iostream>
#include <ctime>
#include "protocol1.h"
#include "protocol2.h"
#include "cpaenc.h"
#include "tools.h"
#include <ctime>

void keygen(unsigned char* pk , unsigned char*sk){
    
    unsigned char seed1[32];
    unsigned char seed2[32];
    unsigned char seed3[32];
    unsigned char k[32];
    randseed(seed1 , 32);
    randseed(seed2 , 32);
    randseed(seed3 , 32);
    randseed(k , 32);
    
    keccak_state state;
    shake256_init(&state);
    shake256_absorb(&state , seed1 , 32);
    polymatq *A = new polymatq(&state ,0, N , N , 1);

    shake256_init(&state);
    shake256_absorb(&state , seed2 , 32);
    polyvecq *b = new polyvecq(&state , 0 , N , 1);

    shake256_init(&state);
    shake256_absorb(&state , seed3 , 32);
    polyvecq *s0 = new polyvecq(&state , ETAR ,N, 0);
    polyvecq *s1 = new polyvecq(&state , ETAR ,N, 0);
    polyvecq *s2 = new polyvecq(&state , ETAR ,N, 0);
    polyvecq *e1 = new polyvecq(&state , ETAR ,N, 0);
    polyvecq *t = new polyvecq(N);

    packsk(sk , s0,s1,s2,e1 , k);

    A->trans();
       
    A->right_mul(t , s2);
    e1->to_ntt();
    t->add(t , e1);
    
    packpk(pk , seed1,seed2, t);

}

void AT_CQ(unsigned char* pi1 ,unsigned char*crsseed, unsigned char* query, unsigned char* H1input ,unsigned char* xchar , unsigned char* pk){

    clock_t t1 = clock();
    unsigned char seed1[32] , tag[MD/8];
    polymatq *A;
    randseed(tag , MD/8);
    memcpy(seed1 , pk , 32);
    keccak_state state;
    shake256_init(&state);
    shake256_absorb(&state , seed1 , 32);
    
    A = new polymatq(&state ,0, N , N , 1);

    shake256_init(&state);
    randseed(seed1 , 32);
    shake256_absorb(&state , seed1 , 32);
    polyvecq *x = new polyvecq(&state , ETAR ,N, 0);
    polyvecq *e2 = new polyvecq(&state , ETAR ,N, 0);

    unsigned char Ginput[VECETALEN*2];
    unsigned char rho[32];
    unsigned char seed2[32];
    
    x->to_char(Ginput , ETAR);
    memcpy(xchar , Ginput , VECETALEN);
    e2->to_char(Ginput+VECETALEN , ETAR);

    shake256_init(&state);
    shake256_absorb(&state , Ginput , VECETALEN*2);
    shake256_squeeze(rho , 32 , &state);

    //unsigned char H1input[256+MD];
    memcpy(H1input , tag , MD/8);
    memcpy(H1input , rho , 32);
    shake256_init(&state);
    shake256_absorb(&state , H1input , MD/8+32);
    polyvecq *ctemp = new polyvecq(&state , 0 ,N, 0);

    ctemp->add(ctemp , e2);
    ctemp->to_ntt();
    polyvecq *c = new polyvecq(N);
    A->right_mul(c , x);
    c->add(c , ctemp);
    delete ctemp;
    c->to_char(query);
    
    //pi1 
    unsigned char pke[64];
    unsigned char cipher[UVECLEN+VPOLYLEN*K3];
    unsigned char uv[POLYPLEN*(K3+4)];
    unsigned char rseed[32];
    unsigned char m[POLYETALEN+32];
    memcpy(m ,xchar , POLYETALEN);
    memcpy(m + POLYETALEN , tag , 32);
    keygen(pke);
    enc(cipher, uv , rseed , pke , m , POLYETALEN+32);
    setup(crsseed);
    x->to_poly();

    clock_t t2 = clock();
    std::cout << "CQ " << (double)(t2-t1)/CLOCKS_PER_SEC << std::endl;

    while(1){
        if(Prove1(pi1 , crsseed , rseed , x , e2 , H1input , uv , pke , A , c) == 1){
            break;
        }
    }
    clock_t t3 = clock();
    std::cout <<"prove1 " << (double)(t3-t2)/CLOCKS_PER_SEC << std::endl;
    std::cout << "done" << std::endl;


}

void AT_IT(unsigned char* pi2 , unsigned char* crsseed2 , unsigned char* resp ,unsigned char* pk, unsigned char *sk , int bit ,\
     unsigned char* query , unsigned char* crsseed , unsigned char* pi1){
    clock_t t1 = clock();
    polyvecq *s0 = new polyvecq(N);
    polyvecq *s1 = new polyvecq(N);
    polyvecq *s2 = new polyvecq(N);
    polyvecq *e1 = new polyvecq(N);
    polyvecq *c = new polyvecq(N);
    unsigned char seed1[32];
    unsigned char seedA[32];
    unsigned char H2input[VECLEN+32];

    memcpy(seedA , pk , 32);
    memcpy(seed1 , pk+32 , 32);
    keccak_state state;
    shake256_init(&state);
    shake256_absorb(&state , seed1 , 32);
    
    polyvecq *b = new polyvecq(&state , 0 , N , 1);

    keccak_state Astate;
    shake256_init(&Astate);
    shake256_absorb(&Astate , seedA , 32);
    
    polymatq *A = new polymatq(&Astate ,0, N , N , 1);

    unpacksk(sk,s0,s1,s2,e1,H2input);
    
    memcpy(H2input+32 , query , VECLEN);

    shake256_init(&state);
    shake256_absorb(&state , H2input , VECLEN+32);

    c->from_char(query , 1);
    Polyq* e3 = new Polyq(&state , ETAE , 0);
    Polyq* h = new Polyq(1);
    Polyq* temppoly = new Polyq(1);

    s2->mul(h , c);
    
    if(bit == 1){
        s1->mul(temppoly , b);
    }
    else{
        s0->mul(temppoly , b);
    }
    h->add(h , temppoly);
    h->to_poly();
    h->add(h,e3);
    h->to_char(resp);
    clock_t t2 = clock();
    std::cout <<"IT " << (double)(t2-t1)/CLOCKS_PER_SEC << std::endl;

    std::cout << "verify1 "<< Verify1(crsseed, pi1 , A , c) << std::endl;
    clock_t t3 = clock();
    std::cout <<"verify1 " << (double)(t3-t2)/CLOCKS_PER_SEC << std::endl;
    setup(crsseed2);
    polyvecq *t = new polyvecq(N);
    t->from_char(pk+64 , 1);
    
    int loop = 0;
    while(!loop){
        if(bit == 1){
            loop = Prove2(pi2 , crsseed2 , s1 , s2 , e1 , e3 , h , A , t , b , c);
        }
        else{
            loop = Prove2(pi2 , crsseed2 , s0 , s2 , e1 , e3 , h , A , t , b , c);
        }
    }
    clock_t t4 = clock();
    std::cout <<"prove2 " << (double)(t4-t3)/CLOCKS_PER_SEC << std::endl;
}

void AT_CF(unsigned char* msg , unsigned char* pk , unsigned char*resp ,unsigned char*H1input, unsigned char* xchar ,\
    unsigned char* query , unsigned char* pi2 , unsigned char* crsseed2){
    
    clock_t t1 = clock();
    polyvecq *t = new polyvecq(N);
    polyvecq *x = new polyvecq(N);
    Polyq* h = new Polyq(0);
    Polyq* v = new Polyq(0);
    Polyq* temp = new Polyq(1);
    unsigned char vchar[COMPOLYLEN];
    
    t->from_char(pk+64 , 1);
    h->from_char(resp , 0);
    x->from_char(xchar , 0 , ETAR);
    t->mul(temp , x);
    temp->to_poly();
    h->sub(v , temp);
    Compress(vchar , v);

    memcpy(msg , H1input ,32+MD/8);
    memcpy(msg+32+MD/8 , vchar , COMPOLYLEN);


    polyvecq *c = new polyvecq(N);
    c->from_char(query , 1);
    keccak_state Astate;
    shake256_init(&Astate);
    shake256_absorb(&Astate , pk , 32);
    keccak_state Bstate;
    shake256_init(&Bstate);
    shake256_absorb(&Bstate , pk+32 , 32);
    
    polyvecq *b = new polyvecq(&Bstate , 0 , N , 1);
    polymatq *A = new polymatq(&Astate ,0, N , N , 1);

    clock_t t2 = clock();
    std::cout <<"CF " << (double)(t2-t1)/CLOCKS_PER_SEC << std::endl;

    std::cout << "verify2 " << Verify2(crsseed2 , pi2 , A , b , c , t , h) << std::endl;
    clock_t t3 = clock();
    std::cout <<"verify2 " << (double)(t3-t2)/CLOCKS_PER_SEC << std::endl;
}
int AT_Rb(unsigned char* pk , unsigned char* sk , unsigned char* msg){

    clock_t t1 = clock();
    keccak_state state;
    int bit = 0;
    unsigned char H1input[32+MD/8];
    unsigned char res[COMPOLYLEN];
    memcpy(H1input , msg , 32+MD/8);
    
    shake256_init(&state);
    shake256_absorb(&state , H1input , 64);
    polyvecq *u = new polyvecq(&state , 0 ,N, 0);

    shake256_init(&state);
    shake256_absorb(&state , pk+32 , 32);
    polyvecq *b = new polyvecq(&state , 0, N ,1);

    Polyq *v_ = new Polyq(0);
    polyvecq *s0 = new polyvecq(N);
    polyvecq *s1 = new polyvecq(N);
    polyvecq *s2 = new polyvecq(N);
    Polyq *u_ = new Polyq(1);
    Polyq *temp = new Polyq(1);
    

    s0->from_char(sk , 0,ETAR);
    s1->from_char(sk+VECETALEN ,0, ETAR);
    s2->from_char(sk+2*VECETALEN,0 , ETAR);
    Decompress(v_ , msg+32+MD/8);
    s2->mul(temp , u);

    s0->mul(u_ , b);
    u_->add(u_ , temp);
    u_->to_poly();
    u_->sub(u_ ,v_);
    
    Compress(res , u_);
    printfBstr(res , COMPOLYLEN);
    for(int i = 0 ; i < COMPOLYLEN;i++){
        if(res[i] != 0){
            bit = 1;
            break;
        }
    }
    if(bit == 0){
        return bit;
    }
    s1->mul(u_ , b);
    u_->add(u_ , temp);
    u_->to_poly();
    u_->sub(u_ ,v_);
    //u_->output();
    Compress(res , u_);
    printfBstr(res , COMPOLYLEN);
    for(int i = 0 ; i < COMPOLYLEN;i++){
        if(res[i] != 0){
            bit = -1;
            break;
        }
    }
    clock_t t2 = clock();
    std::cout <<"Rb " << (double)(t2-t1)/CLOCKS_PER_SEC << std::endl;
    return bit;

}

int main(){
    unsigned char pk[PKLEN] , sk[SKLEN];
    unsigned char query[VECLEN],pi1[PI1LEN] , crsseed[32],H1input[32+MD/8] , xchar[VECETALEN];
    unsigned char resp[VECLEN] , pi2[PI2LEN] , crsseed2[32];
    unsigned char msg[VECLEN+32+MD/8];
    int bit;
    keygen(pk , sk);

    AT_CQ(pi1 ,crsseed, query , H1input , xchar , pk);
    
    AT_IT(pi2 , crsseed2 , resp , pk , sk , 1 , query , crsseed , pi1);
    AT_CF(msg , pk , resp , H1input , xchar , query , pi2 , crsseed2);
    bit = AT_Rb(pk , sk , msg);
    std::cout << bit << std::endl;
}
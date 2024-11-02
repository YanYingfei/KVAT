#include "main_protocol.h"
#include "pack.h"
#include "params.h"

#include <stdlib.h>
#include <iostream>
#include <ctime>

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
    shake256_absorb(&state , seed3 , 32);
    polyvecq *s0 = new polyvecq(&state , ETAR ,N, 0);
    polyvecq *s1 = new polyvecq(&state , ETAR ,N, 0);
    polyvecq *s2 = new polyvecq(&state , ETAR ,N, 0);
    polyvecq *e0 = new polyvecq(&state , ETAR ,N, 0);
    polyvecq *e1 = new polyvecq(&state , ETAR ,N, 0);
    polyvecq *t0 = new polyvecq(N);
    polyvecq *t1 = new polyvecq(N);
    packsk(sk , s0,s1, k);

    A->trans();
       
    A->right_mul(t0 , s0);
    
    e0->to_ntt();
    t0->add(t0 , e0);

    A->right_mul(t1 , s1);
    
    e1->to_ntt();
    t1->add(t1 , e1);
    
    packpk(pk , seed1, t0 , t1);

}

void AT_CQ(unsigned char* pi1 ,unsigned char*crsseed, unsigned char* query, unsigned char* H1input ,unsigned char* xchar , unsigned char* pk, double* pitime){

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
    unsigned char seed2[32];
    
    x->to_char(Ginput , ETAR);
    memcpy(xchar , Ginput , VECETALEN);
    e2->to_char(Ginput+VECETALEN , ETAR);

    //unsigned char H1input[256+MD];
    memcpy(H1input , tag , MD/8);
    shake256_init(&state);
    shake256_absorb(&state , H1input , MD/8);
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
    // std::cout << "CQ " << (double)(t2-t1)/CLOCKS_PER_SEC << std::endl;
    int cnt=0;
    while(1){
        cnt++;
        if(Prove1(pi1 , crsseed , rseed , x , e2 , H1input , uv , pke , A , c) == 1){
            std::cout<< "loop time:"<< cnt<< std::endl;
            break;
        }
    }
    clock_t t3 = clock();
    if(pitime != NULL){
        pitime[0] = (double)(t3-t2)/CLOCKS_PER_SEC;
    }
    // std::cout <<"prove1 " << (double)(t3-t2)/CLOCKS_PER_SEC << std::endl;
    std::cout << "done" << std::endl;
}

void AT_IT(unsigned char* pi2 , unsigned char* crsseed2 , unsigned char* resp ,unsigned char* pk, unsigned char *sk ,\
     unsigned char* query , unsigned char* crsseed , unsigned char* pi1, double* pi1time , double *pi2time){
    clock_t time1 = clock();
    polyvecq *s0 = new polyvecq(N);
    polyvecq *s1 = new polyvecq(N);
    polyvecq *s2 = new polyvecq(N);
    polyvecq *eb = new polyvecq(N);
    polyvecq *c = new polyvecq(N);
    polyvecq *t0 = new polyvecq(N);
    polyvecq *t1 = new polyvecq(N);
    int bit;
    unsigned char seed1[32];
    unsigned char seedA[32];
    unsigned char H2input[VECLEN+32];

    memcpy(seedA , pk , 32);
    memcpy(seed1 , pk+32 , 32);
    keccak_state state;

    keccak_state Astate;
    shake256_init(&Astate);
    shake256_absorb(&Astate , seedA , 32);
    
    polymatq *A = new polymatq(&Astate ,0, N , N , 1);

    t0->from_char(pk+32 , 1);
    t1->from_char(pk+32+VECLEN , 1);
    unpacksk(sk,s0,s1,H2input);
    
    memcpy(H2input+32 , query , VECLEN);

    shake256_init(&state);
    shake256_absorb(&state , H2input , VECLEN+32);

    c->from_char(query , 1);
    Polyq* e3 = new Polyq(&state , ETAE , 0);
    Polyq* h = new Polyq(1);
    Polyq* temppoly = new Polyq(1);

    bit = rand() % 2;
    if(bit == 1){
        s1->mul(h , c);
        A->trans();
        A->right_mul(eb , s1);
        t1->sub(eb , eb);
        eb->to_poly();
    }
    else{
        s0->mul(h , c);
        A->trans();
        A->right_mul(eb , s0);
        t0->sub(eb , eb);
        eb->to_poly();
    }
    A->trans();
    h->to_poly();
    h->add(h,e3);
    h->to_char(resp);
    clock_t t2 = clock();
    std::cout << "Private bit in TokenIssue: " << bit << std::endl;
    // std::cout <<"IT " << (double)(t2-time1)/CLOCKS_PER_SEC << std::endl;

    std::cout << "verify Pi1: "<< Verify1(crsseed, pi1 , A , c) << std::endl;
    clock_t t3 = clock();
    if(pi1time != NULL){
        pi1time[0] = (double)(t3-t2)/CLOCKS_PER_SEC;
    }
    // std::cout <<"verify1 " << (double)(t3-t2)/CLOCKS_PER_SEC << std::endl;
    setup(crsseed2);
    // polyvecq *t = new polyvecq(N);
    // t->from_char(pk+64 , 1);
    
    int loop = 0;
    while(!loop){
        if(bit == 1){
            loop = Prove2(pi2 , crsseed2 , 1 , s1 , eb , e3 , h , A , t0, t1 , c);
        }
        else{
            loop = Prove2(pi2 , crsseed2 , 0 , s0 , eb , e3 , h , A , t0 ,t1 , c);
        }
        // std::cout << "1"<<std::endl;
    }
    clock_t t4 = clock();
    // if(pi2time != NULL){
    //     pi2time[0] = (double)(t4-t3)/CLOCKS_PER_SEC;
    // }
    // std::cout <<"prove2 " << (double)(t4-t3)/CLOCKS_PER_SEC << std::endl;
}

void AT_CF(unsigned char* msg , unsigned char* pk , unsigned char*resp ,unsigned char*H1input, unsigned char* xchar ,\
    unsigned char* query , unsigned char* pi2 , unsigned char* crsseed2, double* pitime){
    
    clock_t time1 = clock();
    polyvecq *t0 = new polyvecq(N);
    polyvecq *t1 = new polyvecq(N);
    polyvecq *x = new polyvecq(N);
    Polyq* h = new Polyq(0);
    Polyq* v0 = new Polyq(0);
    Polyq* v1 = new Polyq(0);
    Polyq* temp = new Polyq(1);
    unsigned char v0char[COMPOLYLEN] , v1char[COMPOLYLEN];
    
    t0->from_char(pk+32 , 1);
    t1->from_char(pk+32+VECLEN , 1);
    h->from_char(resp , 0);
    x->from_char(xchar , 0 , ETAR);
    t0->mul(temp , x);
    temp->to_poly();
    h->sub(v0 , temp);
    Compress(v0char , v0);

    t1->mul(temp , x);
    temp->to_poly();
    h->sub(v1 , temp);
    Compress(v1char , v1);

    memcpy(msg , H1input ,32);
    memcpy(msg+32 , v0char , COMPOLYLEN);
    memcpy(msg+32+COMPOLYLEN , v1char , COMPOLYLEN);


    polyvecq *c = new polyvecq(N);
    c->from_char(query , 1);
    keccak_state Astate;
    shake256_init(&Astate);
    shake256_absorb(&Astate , pk , 32);

    polymatq *A = new polymatq(&Astate ,0, N , N , 1);

    clock_t t2 = clock();
    // std::cout <<"CF " << (double)(t2-time1)/CLOCKS_PER_SEC << std::endl;

    std::cout << "verify Pi2: " << Verify2(crsseed2 , pi2 , A , c , t0 , t1 , h) << std::endl;
    clock_t t3 = clock();
    if(pitime != NULL){
        pitime[0] = (double)(t3-t2)/CLOCKS_PER_SEC;
    }
    // std::cout <<"verify2 " << (double)(t3-t2)/CLOCKS_PER_SEC << std::endl;
}
int AT_Rb(unsigned char* pk , unsigned char* sk , unsigned char* msg){
    std::cout <<"Read bit: " << std::endl;
    clock_t t1 = clock();
    keccak_state state;
    int bit = 0;
    unsigned char H1input[32];
    unsigned char res[COMPOLYLEN];
    memcpy(H1input , msg , 32);
    
    shake256_init(&state);
    shake256_absorb(&state , H1input , 32);
    polyvecq *h = new polyvecq(&state , 0 ,N, 0);
    Polyq *v0_ = new Polyq(0);
    Polyq *v1_ = new Polyq(0);
    polyvecq *s0 = new polyvecq(N);
    polyvecq *s1 = new polyvecq(N);
    Polyq *u0 = new Polyq(1);
    Polyq *u1 = new Polyq(1);
    Polyq *temp = new Polyq(1);
    

    s0->from_char(sk , 0,ETAR);
    s1->from_char(sk+VECETALEN ,0, ETAR);
    Decompress(v0_ , msg+32);
    Decompress(v1_ , msg+32+COMPOLYLEN);
    
    s0->mul(u0 , h);
    s1->mul(u1 , h);
    
    u0->to_poly();
    u1->to_poly();
    
    u0->sub(temp , v0_);
    Compress(res , temp);
    //printfBstr(res , COMPOLYLEN);
    for(int i = 0 ; i < COMPOLYLEN;i++){
        if(res[i] != 0){
            bit = 1;
            break;
        }
    }
    if(bit == 0){
        return bit;
    }
    u1->sub(temp ,v1_);
    //u_->output();
    Compress(res , temp);
    // printfBstr(res , COMPOLYLEN);
    for(int i = 0 ; i < COMPOLYLEN;i++){
        if(res[i] != 0){
            bit = -1;
            break;
        }
    }
    clock_t t2 = clock();
    // std::cout <<"Rb " << (double)(t2-t1)/CLOCKS_PER_SEC << std::endl;

    return bit;

}


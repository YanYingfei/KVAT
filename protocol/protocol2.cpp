#include <stdlib.h>
#include <iostream>
#include <ctime>
#include "../poly/poly.h"
#include "../poly/polymat.h"
#include "../poly/sample.h"
#include "params.h"
#include "tools.h"
#include "construct.h"
#include <cmath>
#include "protocol1.h"

void generatem2(Polyq *m1list[] ,polyvecq* sb ,polyvecq* s2,polyvecq* e1,Polyq *e3){
    
    sb->to_poly();
    s2->to_poly();
    e1->to_poly();
    for(int i = 0;i<N;i++){
        m1list[i] = new Polyq(0);
        m1list[i]->copy(sb->polyarray[i]);
        m1list[i+1+N] = new Polyq(0);
        m1list[N+1+i]->copy(s2->polyarray[i]);
        m1list[N*2+2+i] = new Polyq(0);
        m1list[N*2+2+i]->copy(e1->polyarray[i]);
    }
    
    m1list[N] = new Polyq(0);
    geta(m1list[N] , BETAE - sb->norm());
    m1list[N*2+1] = new Polyq(0);
    geta(m1list[N*2+1] ,BETAE - s2->norm());
    m1list[N*3+2] = new Polyq(0);
    geta(m1list[N*3+2] ,BETAE - e1->norm());
    //tag 
    m1list[N*3+3] = new Polyq(0);
    m1list[N*3+3]->copy(e3);
    m1list[N*3+4] = new Polyq(0);
    geta(m1list[N*3+4] , BETAE2-e3->norm());

}
void constructFy2(polymatq *Fymat ,polymatq *Fg ,polymatq *F3){

    for(int i = 0 ; i < TAU ; i++){
        Fymat->vecarray[i] = Fg->vecarray[i];
    }
    for(int i = 0 ; i < 2;i++){
        Fymat->vecarray[TAU+i] = F3->vecarray[i];
    }
}
void constructy2(polyvecq *y ,polymatq *Fymat ,polyvecq *y1 ,polyvecq*y2){
    polyvecq temp(y1->k);
    y1->sigma(&temp);
    for(int i = 0 ; i < y1->k ; i++){
        y->polyarray[i]->copy(y1->polyarray[i]);
        y->polyarray[i+y1->k]->copy(temp.polyarray[i]);
    }
    polyvecq temp1(Fymat->k);
    Fymat->right_mul(&temp1 , y2);
    polyvecq temp2(Fymat->k);
    temp1.sigma(&temp2);
    for(int i = 0 ; i < temp1.k ; i++){
        temp1.polyarray[i]->mul_num(y->polyarray[i+y1->k*2] , -1);
        temp2.polyarray[i]->mul_num(y->polyarray[i+y1->k*2+temp1.k] , -1);
    }
}
void computeh(polyvecq* hvec ,polyvecq*g,Polyq *F[],polymatq *Gamma){
    Polyq temp(1);
    hvec->to_ntt();
    for(int j = 0 ; j < TAU ; j++){
        for(int i = 0 ; i < 4 ; i++){
            Gamma->vecarray[j]->polyarray[N+1+i]->mul(&temp , F[i]);
            temp.add(hvec->polyarray[j] , hvec->polyarray[j]);
        }
    }
    hvec->to_poly();
    //hvec->add(hvec , g);
}
void constructD2_2(Polyq *D2list[] ,polyvecq *mu ,polymatq *Gamma){
    
    Polyq temp(1);
    for(int i = 0 ; i < 4;i++){
        D2list[i] = new Polyq(1);
        for(int j = 0 ; j < TAU ; j++){
            mu->polyarray[j]->mul(&temp , Gamma->vecarray[j]->polyarray[N+1+i]);
            D2list[i]->add(D2list[i] , &temp);
        }
    }
}
void constructd1_2(polyvecq *d1vec ,polyvecq* mu ,polymatq *Gamma ,polymatq *A , polyvecq *b , polyvecq *c){
    
    Polyq temp1(1) , temp2(1);
    polyvecq temp3(N);
    for(int j = 0 ; j < TAU ; j++){
        mu->polyarray[j]->mul(&temp1 , Gamma->vecarray[j]->polyarray[0]);
        temp2.add(&temp2 , &temp1);
    }
    d1vec->polyarray[3*N+3]->copy(&temp2);
    b->mul_poly(&temp3 , &temp2);
    for(int i = 0 ; i < N ; i++){
        d1vec->polyarray[i]->copy(temp3.polyarray[i]);
    }

    polyvecq temp4(N) , temp5(N , 1);
    A->trans();
    for(int j = 0 ; j < TAU ; j++){
        temp4.reset(1);
        for(int i = 0 ; i < N ; i++){
            A->vecarray[i]->mul_poly(&temp3 , Gamma->vecarray[j]->polyarray[i+1]);
            temp4.add(&temp4 , &temp3);
        }
        c->mul_poly(&temp3 , Gamma->vecarray[j]->polyarray[0]);
        temp4.add(&temp4 , &temp3);
        
        temp4.mul_poly(&temp4 , mu->polyarray[j]);
        temp5.add(&temp5 , &temp4);
    }
    A->trans();
    for(int i = 0 ; i < N ; i++){
        d1vec->polyarray[i+N+1]->copy(temp5.polyarray[i]);
    }
    temp5.reset(1);
    for(int j = 0 ; j < TAU ; j++){
        
        for(int i = 0 ; i < N;i++){
            temp4.polyarray[i]->copy(Gamma->vecarray[j]->polyarray[i+1]);
        }
        temp4.mul_poly(&temp4 , mu->polyarray[j]);
        temp5.add(&temp5 , &temp4);
    }
    for(int i = 0 ; i < N ; i++){
        d1vec->polyarray[i+2*N+2]->copy(temp5.polyarray[i]);
    }

}
void constructm2(polyvecq *m ,polyvecq *m1 ,polyvecq *g ,polyvecq *y3){
    int mlen = 3*N+5;
    int m2len = TAU+2;
    for(int i = 0 ; i < mlen;i++){
        m->polyarray[i]->copy(m1->polyarray[i]);
        m1->polyarray[i]->sigma(m->polyarray[i+mlen]);
    }
    for(int i = 0 ; i < TAU; i++){
        m->polyarray[2*mlen+i]->copy(g->polyarray[i]);
        g->polyarray[i]->sigma(m->polyarray[i+mlen*2+m2len]);
    }
    for(int i = 0 ; i < 2; i++){
        m->polyarray[2*mlen+TAU+i]->copy(y3->polyarray[i]);
        y3->polyarray[i]->sigma(m->polyarray[i+mlen*2+m2len+TAU]);
    }
}
void constructz2(polyvecq *z ,polymatq*Fymat, polyvecq *tg, polyvecq *t3,Polyq*cpoly, polyvecq *z1, polyvecq *z2){
    polyvecq ty(TAU+2);
    for(int i = 0 ; i < TAU; i++){
        ty.polyarray[i]->copy(tg->polyarray[i]);
    }
    for(int i = 0 ; i < 2; i++){
        ty.polyarray[i+TAU]->copy(t3->polyarray[i]);
    }
    polyvecq temp(z1->k);
    z1->sigma(&temp);
    for(int i = 0 ; i < z1->k ; i++){
        z->polyarray[i]->copy(z1->polyarray[i]);
        z->polyarray[i+z1->k]->copy(temp.polyarray[i]);
    }
    polyvecq temp1(Fymat->k);
    Fymat->right_mul(&temp1 , z2);
    polyvecq temp2(Fymat->k);
    ty.mul_poly(&ty , cpoly);
    ty.sub(&temp1 , &temp1);
    temp1.sigma(&temp2);
    for(int i = 0 ; i < temp1.k ; i++){
        z->polyarray[i+z1->k*2]->copy(temp1.polyarray[i]);
        z->polyarray[i+z1->k*2+temp1.k]->copy(temp2.polyarray[i]);
    }
}
void D2mul2(Polyq *res , Polyq* D2list[] , polyvecq *left , polyvecq *right){
    Polyq temp(1);
    int Len = 3*N+5;
    int Len1 = N+1;
    res->reset(1);
    for(int i = 0 ; i < Len1 ; i++){
        left->polyarray[Len+i]->mul(&temp , D2list[0]);
        right->polyarray[i]->mul(&temp , &temp);
        res->add(res , &temp);
    }
    for(int i = Len1 ; i < 2*Len1 ; i++){
        left->polyarray[Len+i]->mul(&temp , D2list[1]);
        right->polyarray[i]->mul(&temp , &temp);
        res->add(res , &temp);
    }
    for(int i = 2*Len1 ; i < 3*Len1 ; i++){
        left->polyarray[Len+i]->mul(&temp , D2list[2]);
        right->polyarray[i]->mul(&temp , &temp);
        res->add(res , &temp);
    }
    for(int i = 3*Len1 ; i < 3*Len1+2 ; i++){
        left->polyarray[Len+i]->mul(&temp , D2list[3]);
        right->polyarray[i]->mul(&temp , &temp);
        res->add(res , &temp);
    }
}
void constructd2(Polyq *d, polyvecq *mu, polymatq *Gamma, polyvecq *t , Polyq*h , polyvecq* hvec){
    Polyq temp1(1) , temp2(1);
    Polyq temp3(1) , temp4(1);
    d->to_ntt();
    
    hvec->mul(d , mu);
    for(int j = 0 ; j < TAU ; j++){
        mu->polyarray[j]->mul(&temp1 , Gamma->vecarray[j]->polyarray[0]);
        temp2.add(&temp2 , &temp1);
    }
    temp2.mul(&temp2 , h);
    d->add(d , &temp2);
    

    

    for(int j = 0 ; j < TAU; j++){
        temp1.reset(1);
        for(int i = 0 ; i < N ;i++){
            Gamma->vecarray[j]->polyarray[i+1]->mul(&temp2 , t->polyarray[i]);
            temp1.add(&temp1 , &temp2);
        }
        temp1.mul(&temp1 , mu->polyarray[j]);
        d->add(d , &temp1);
    }
    
    temp2.reset(1);
    for(int j = 0 ; j < TAU; j++){
        temp1.reset(1);
        for(int i = 0 ; i < 3 ;i++){
            temp1.add(&temp1 , Gamma->vecarray[j]->polyarray[i+N+1]);
        }
        temp1.mul(&temp1 , mu->polyarray[j]);
        temp2.add(&temp2 , &temp1);
    }
    temp2.mul_num(&temp2 , BETAE);
    d->to_poly();
    d->add(d , &temp2);
    
    for(int j = 0 ; j < TAU; j++){
        Gamma->vecarray[j]->polyarray[N+4]->mul(&temp1 , mu->polyarray[j]);
        temp1.mul_num(&temp2 , BETAE2);
        
        d->add(d , &temp2);
    }

    d->mul_num(d , -1);
}


int Prove2(unsigned char* pi , unsigned char* crsseed ,polyvecq *sb , polyvecq *s2 , polyvecq *e1 , Polyq *e3, Polyq* h,polymatq*A,polyvecq*t, polyvecq *b,\
    polyvecq *c){
    //unpack crs
    keccak_state crsstate;
    shake256_init(&crsstate);
    randseed(crsseed , 32);
    shake256_absorb(&crsstate , crsseed,32);

    polymatq *E1 = new polymatq(&crsstate, 0, N1 ,3*N+5 , 1);
    polymatq *E2 = new polymatq(&crsstate, 0, N1 ,M2 , 1);
    polymatq *Fg = new polymatq(&crsstate, 0, TAU ,M2, 1);
    polymatq *F3 = new polymatq(&crsstate, 0, 2 ,M2 , 1);
    polymatq *F4 = new polymatq(&crsstate, 0, 2 ,M2 , 1);
    polyvecq *fg = new polyvecq(&crsstate, 0, M2, 1);
    // generate m
    Polyq *m1list[3*N+5];
    generatem2(m1list , sb , s2,e1,e3);
    
    polyvecq *m1 = new polyvecq(m1list , 3*N+5);

    //generate e

    unsigned char random_seed[32];
    keccak_state random_state;
    shake256_init(&random_state);
    randseed(random_seed , 32);
    shake256_absorb(&random_state , random_seed,32);
    
    polyvecq *e = new polyvecq(&random_state , ETAR , M2 , 0);


    // FR step 1

    polyvecq temp1(N1);
    polyvecq temp2(N1);
    polyvecq *tE = new polyvecq(N1);
    
    E1->right_mul(&temp1 ,m1);
    E2->right_mul(&temp2 , e);
    
    temp1.add(tE , &temp2);
    tE->to_poly();
    //FR step 2
    polyvecq *y1 = new polyvecq(Y1 , 3*N+5,0);
    polyvecq *y2 = new polyvecq(Y2 , M2,0);
    polyvecq *w = new polyvecq(N1);

    E1->right_mul(&temp1 , y1);
    E2->right_mul(&temp2 , y2);
    temp1.add(w , &temp2);
    w->to_poly();

    //FR step 3
    polyvecq *g = new polyvecq(&random_state,0 ,TAU , 0);
    polyvecq *tg = new polyvecq(TAU);
    for(int i = 0;i< TAU;i++){
        g->polyarray[i]->polyarray[0] = 0;
    }
    Fg->right_mul(tg , e);
    tg->to_poly();
    tg->add(tg , g);

    //FR step 4-5
    polyvecq *y3 = new polyvecq(Y3 , 2 , 0);
    
    polyvecq *t3 = new polyvecq(2);
    polyvecq *t4 = new polyvecq(2);
    
    F3->right_mul(t3 , e);
    t3->to_poly();
    t3->add(t3 , y3);
    
    //FR step 6
    F4->right_mul(t4 , e);
    t4->to_poly();

    int bits[2];
    for(int i = 0 ; i < 2;i++){
        bits[i] = 2*(rand()%2)-1;
        t4->polyarray[i]->polyarray[0] =t4->polyarray[0]->mod(t4->polyarray[i]->polyarray[0]+bits[i]);
    }

    //SR
    unsigned char Rseed[POLYLEN*(4+TAU+2*N1)+32+VECETALEN];//VECETALEN = XSTATEMENTLEN
    memcpy(Rseed , crsseed , 32);

    // todo: insert statement x
    memset(Rseed+32 , 0, VECETALEN);

    t3->to_char(Rseed+32+VECETALEN , 0);
    t4->to_char(Rseed+32+VECETALEN+POLYLEN*2 , 0);
    tg->to_char(Rseed+32+VECETALEN+POLYLEN*4 , 0);
    tE->to_char(Rseed+32+VECETALEN+POLYLEN*(4+TAU) , 0);
    w->to_char(Rseed+32+VECETALEN+POLYLEN*(4+TAU+N1) , 0);

    keccak_state Rstate;
    shake256_init(&Rstate);
    shake256_absorb(&Rstate , Rseed,POLYLEN*(4+TAU+2*N1)+32+VECETALEN);
    
    int64_t *R = (int64_t*)malloc(sizeof(int64_t)*256*(3*N+5)*128);

    cbd(&Rstate , R , 1 , 256*128*(3*N+5));

    //SR step4
    int64_t *m1vec = (int64_t*)malloc(sizeof(int64_t)*(3*N+5)*128);
    int64_t z3vec[256];
    int64_t Rm1[256];
    int64_t y3vec[256];
    polyvec2vec(m1vec , m1);
    polyvec2vec(y3vec , y3);
    matrix_mul(Rm1 , R , m1vec , 256 , (3*N+5)*128);
    memcpy(z3vec , Rm1 , 256*sizeof(int64_t));
    vec_mul_si(z3vec , bits[0] , 256);
    vec_add(z3vec , y3vec , 256);
    
    if(Rej0(z3vec ,Rm1 , Y3 , 7)){
        return 0;
    }


    // TR
    unsigned char Gammaseed[256*8];
    vec2char(Gammaseed , z3vec , 256);

    keccak_state gamma_state;
    shake256_init(&gamma_state);
    shake256_absorb(&gamma_state , Rseed,POLYLEN*(4+TAU+2*N1)+32+VECETALEN);
    shake256_absorb(&gamma_state , Gammaseed , 256*8);
    
    polymatq *Gamma = new polymatq(TAU , N+5);
    sample_gamma(&gamma_state , Gamma);

    // compute F
    Polyq *F[4];


    polyvecq *temp8 = new polyvecq(N+1);
    polyvecq *temp9 = new polyvecq(N+1);
    int tempbetalist[4] = {BETAE , BETAE , BETAE , BETAE2};
    for(int j = 0 ; j < 3 ; j++){
        F[j] = new Polyq(0);
        for(int i = 0 ; i < N+1 ; i++){
            temp8->polyarray[i]->copy(m1list[i+(N+1)*j]);
        }
        temp8->sigma(temp9);
        temp8->mul(F[j] , temp9);
        F[j]->to_poly();
        F[j]->polyarray[0] -= tempbetalist[j];
    }
    polyvecq *temp8_1 = new polyvecq(2);
    polyvecq *temp9_1 = new polyvecq(2);
    for(int i = 0 ; i < 2 ; i++){
        temp8_1->polyarray[i]->copy(m1list[3*(N+1)+i]);
    }
    F[3] = new Polyq(0);
    temp8_1->sigma(temp9_1);
    temp8_1->mul(F[3] , temp9_1);
    F[3]->to_poly();
    F[3]->polyarray[0] -= tempbetalist[3];
    

    //compute hvec
    polyvecq *hvec = new polyvecq(TAU);
    computeh(hvec , g,F,Gamma);
    //FourR step1

    //H3input

    unsigned char museed[TAU*POLYLEN];
    hvec->to_char(museed);
    keccak_state mu_state;
    shake256_init(&mu_state);
    shake256_absorb(&mu_state , Rseed,POLYLEN*(4+TAU+2*N1)+32+VECETALEN);
    shake256_absorb(&mu_state , Gammaseed , 256*8);
    shake256_absorb(&mu_state , museed , TAU*POLYLEN);

    polyvecq *mu = new polyvecq(&mu_state ,0, TAU , 0);


    //FourR step 2
    polymatq *Fymat = new polymatq(TAU+2 , M2);
    polyvecq *y = new polyvecq((3*N+5)*2 + (TAU+2)*2);
    constructFy2(Fymat ,Fg , F3);
    constructy2(y , Fymat , y1 , y2);

    //FourR step 3
    Polyq *D2list[4];
    constructD2_2(D2list , mu , Gamma);

    polyvecq *d1vec = new polyvecq((3*N+5)*2 + (TAU+2)*2);
    constructd1_2(d1vec , mu , Gamma , A , b , c);

    //FourR step 4
    
    //construct m
    polyvecq *m = new polyvecq((3*N+5)*2 + (TAU+2)*2);
    constructm2(m , m1 , g , y3 );

    //compute g1
    Polyq tempg(1);
    Polyq *g1 = new Polyq(1);
    D2mul2(&tempg , D2list , m , y);
    g1->add(g1 , &tempg);
    D2mul2(&tempg , D2list,y , m);
    g1->add(g1 , &tempg);
    d1vec->mul(&tempg , y);
    g1->add(g1 , &tempg);

    
    //compute t
    Polyq *tpoly = new Polyq(1);
    fg->mul(tpoly , e);
    tpoly->add(tpoly , g1);
    tpoly->to_poly();

    //compute v
    Polyq *v = new Polyq(1);
    D2mul2(&tempg , D2list , y ,y);
    v->add(v , &tempg);
    fg->mul(&tempg , y2);
    v->add(v , &tempg);
    v->to_poly();

    //FifthR

    //H4input
    unsigned char cseed[POLYLEN*2];

    tpoly->to_char(cseed);
    v->to_char(cseed+POLYLEN);
    keccak_state c_state;
    shake256_init(&c_state);
    shake256_absorb(&c_state , Rseed,POLYLEN*(4+TAU+2*N1)+32+VECETALEN);
    shake256_absorb(&c_state , Gammaseed , 256*8);
    shake256_absorb(&c_state , museed , TAU*POLYLEN);
    shake256_absorb(&c_state , cseed , POLYLEN*2);
    Polyq *cpoly = new Polyq(0);
    Csamplepoly(&c_state , cpoly);
    polyvecq *cm1 = new polyvecq(m1->k);
    polyvecq *z1 = new polyvecq(m1->k);
    m1->mul_poly(cm1 , cpoly);
    cm1->to_poly();
    z1->add(z1 , cm1);
    y1->to_poly();
    z1->add(z1 , y1);

    polyvecq *ce = new polyvecq(e->k);
    polyvecq *z2 = new polyvecq(e->k);
    //cpoly->mul_num(cpoly , bits[1]);
    e->mul_poly(ce , cpoly);
    ce->to_poly();
    z2->add(z2 , ce);
    y2->to_poly();
    z2->add(z2 , y2);

    if(Rej1(z1 , cm1 , Y1 , 18) == 0){
        return 0;
        //goto reboot;
    }
    if(Rej1(z2,ce , Y2 , 15) == 0){
        return 0;
        //goto reboot;
    }
    
    memcpy(pi , Rseed+VECETALEN+32 , POLYLEN*(4+TAU+2*N1));//POLYLEN*(7+TAU+N1*2+N)
    memcpy(pi+POLYLEN*(4+TAU+2*N1) , Gammaseed,256*8);
    memcpy(pi + POLYLEN*(4+TAU+2*N1) + 256*8 , museed, TAU*POLYLEN);

    
    memcpy(pi + POLYLEN*(4+TAU*2+2*N1) + 256*8 , cseed , POLYLEN*2);
//2*N+M+K1+3
    z1->to_char(pi + POLYLEN*(6+TAU*2+2*N1) + 256*8);
    z2->to_char(pi + POLYLEN*(11+TAU*2+2*N1+3*N) + 256*8);


    // 16+TAU*2+2*N1+6*N
    return 1;
}

int Verify2(unsigned char *crsseed, unsigned char *pi , polymatq *A, polyvecq *b, polyvecq *c , polyvecq *t , Polyq *h){
    keccak_state crsstate;
    shake256_init(&crsstate);
    shake256_absorb(&crsstate , crsseed,32);

    polymatq *E1 = new polymatq(&crsstate, 0, N1 ,3*N+5 , 1);
    polymatq *E2 = new polymatq(&crsstate, 0, N1 ,M2 , 1);
    polymatq *Fg = new polymatq(&crsstate, 0, TAU ,M2, 1);
    polymatq *F3 = new polymatq(&crsstate, 0, 2 ,M2 , 1);
    polymatq *F4 = new polymatq(&crsstate, 0, 2 ,M2 , 1);
    polyvecq *fg = new polyvecq(&crsstate, 0, M2, 1);
    
    unsigned char Rseed[POLYLEN*(7+TAU+N1*2+N)+32+VECETALEN];
    memcpy(Rseed , crsseed , 32);

    //todo: insert statement x
    memset(Rseed+32 , 0 , VECETALEN);

    memcpy(Rseed+VECETALEN+32 , pi , POLYLEN*(7+TAU+N1*2+N));

    // unpack alpha1
    polyvecq *t3 = new polyvecq(2);
    polyvecq *t4 = new polyvecq(2);
    polyvecq *tg = new polyvecq(TAU);
    polyvecq *tE = new polyvecq(N1);
    polyvecq *w = new polyvecq(N1);
    t3->from_char(Rseed+32+VECETALEN , 0);
    t4->from_char(Rseed+32+VECETALEN+POLYLEN*2 , 0);
    tg->from_char(Rseed+32+VECETALEN+POLYLEN*4 , 0);
    tE->from_char(Rseed+32+VECETALEN+POLYLEN*(4+TAU) , 0);
    w->from_char(Rseed+32+VECETALEN+POLYLEN*(4+TAU+N1) , 0);


    keccak_state Rstate;
    shake256_init(&Rstate);
    shake256_absorb(&Rstate , Rseed,POLYLEN*(4+TAU+2*N1)+32+VECETALEN);
    
    int64_t *R = (int64_t*)malloc(sizeof(int64_t)*256*(3*N+5)*128);

    cbd(&Rstate , R , 1 , 256*128*(3*N+5));
    clock_t t1 = clock();
    
    //unpack alpha2

    int64_t z3vec[256];
    unsigned char Gammaseed[256*8];
    memcpy(Gammaseed ,pi+POLYLEN*(4+TAU+2*N1) ,256*8);
    char2vec(z3vec , Gammaseed , 256);

    keccak_state gamma_state;
    shake256_init(&gamma_state);
    shake256_absorb(&gamma_state , Rseed,POLYLEN*(4+TAU+2*N1)+32+VECETALEN);
    shake256_absorb(&gamma_state , Gammaseed , 256*8);
    
    polymatq *Gamma = new polymatq(TAU , N+5);
    sample_gamma(&gamma_state , Gamma);


    // unpack alpha3
    unsigned char museed[TAU*POLYLEN];
    memcpy(museed, pi + POLYLEN*(4+TAU+2*N1) + 256*8 , TAU*POLYLEN);
    polyvecq *hvec = new polyvecq(TAU);
    hvec->from_char(museed , 0);
    keccak_state mu_state;
    shake256_init(&mu_state);
    shake256_absorb(&mu_state , Rseed,POLYLEN*(4+TAU+2*N1)+32+VECETALEN);
    shake256_absorb(&mu_state , Gammaseed , 256*8);
    shake256_absorb(&mu_state , museed , TAU*POLYLEN);

    polyvecq *mu = new polyvecq(&mu_state ,0, TAU , 0);

    //unpack alpha4
    unsigned char cseed[POLYLEN*2];
    memcpy(cseed , pi + POLYLEN*(4+TAU*2+2*N1) + 256*8 , POLYLEN*2);
    Polyq *tpoly = new Polyq(0);
    Polyq *v = new Polyq(0);
    tpoly->from_char(cseed , 0);
    v->from_char(cseed+POLYLEN , 0);
    keccak_state c_state;
    shake256_init(&c_state);
    shake256_absorb(&c_state , Rseed,POLYLEN*(4+TAU+2*N1)+32+VECETALEN);
    shake256_absorb(&c_state , Gammaseed , 256*8);
    shake256_absorb(&c_state , museed , TAU*POLYLEN);
    shake256_absorb(&c_state , cseed , POLYLEN*2);
    Polyq *cpoly = new Polyq(0);
    Csamplepoly(&c_state , cpoly);
    //constructz
    polyvecq *z1 = new polyvecq(3*N+5);
    polyvecq *z2 = new polyvecq(M2);
    polyvecq *z = new polyvecq(2*(3*N+5)+2*(2+TAU));
    z1->from_char(pi + POLYLEN*(6+TAU*2+2*N1) + 256*8 , 0);
    z2->from_char(pi + POLYLEN*(11+TAU*2+2*N1+3*N) + 256*8 , 0);
    
    polymatq *Fymat = new polymatq(TAU+2 , M2);
    constructFy2(Fymat ,Fg , F3);
    constructz2(z , Fymat , tg , t3,cpoly , z1 , z2);
    // constructD
    Polyq *D2list[4];
    constructD2_2(D2list , mu , Gamma);
    polyvecq *d1vec = new polyvecq(2*(3*N+5)+2*(2+TAU));
    constructd1_2(d1vec , mu , Gamma , A , b , c);

    Polyq *d0 = new Polyq(1);
    constructd2(d0 , mu , Gamma , t , h , hvec);
    //check norm
    if(z1->norm() > B1*B1 || z2->norm() > B2*B2 || vec_norm(z3vec , 256) > B3*B3){
        std::cout << "wrong1" << std::endl;
        return 0;
    }

    // check E1z1 + E2z2 == w+ ctE

    polyvecq left(N1 , 1) , temp(N1 , 1) , right(N1 , 1);
    E1->right_mul(&temp , z1);
    left.add(&left, &temp);
    E2->right_mul(&temp , z2);
    left.add(&left , &temp);

    w->to_ntt();
    tE->mul_poly(&right , cpoly);
    right.add(&right , w);
    

    if(left.equal(&right) == 0){
        std::cout << "wrong2" << std::endl;
        return 0;
    }
    // check z^TD2z + cd1^Tz + c^2d0 + fg^Tz2= v + ct

    Polyq leftpoly(1) , temppoly(1);
    
    D2mul2(&leftpoly , D2list ,z ,  z);
    
    d1vec->mul(&temppoly , z);
    temppoly.mul(&temppoly , cpoly);
    leftpoly.add(&leftpoly , &temppoly);
    
    d0->mul(&temppoly , cpoly);
    
    cpoly->mul(&temppoly , &temppoly);

    leftpoly.add(&leftpoly , &temppoly);
    fg->mul(&temppoly , z2);
    leftpoly.add(&leftpoly , &temppoly);

    cpoly->mul(&temppoly , tpoly);
    leftpoly.sub(&leftpoly , &temppoly); 
    if(leftpoly.equal(v) == 0){
        std::cout << "wrong3" << std::endl;
        return 0;
    }

    //check h
    hvec->to_poly();
    for(int i = 0 ; i < TAU ;i++){
        if(hvec->polyarray[i]->polyarray[0] !=0 ){
            std::cout << "wrong4" << std::endl;
            return 0;
        }
    }
   
    return 1;
}
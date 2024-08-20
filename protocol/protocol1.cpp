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
void setup(unsigned char* crsseed){
    randseed(crsseed , 32);
}



void encodem(Polyq* mpoly , unsigned char*m){
    for(int i = 0 ; i < 16;i++){
        for(int j = 7 ; j >=0;j--){
            if((m[i]>>j)&1 == 1){
                mpoly->polyarray[i*8+j] += 3329/2;
            }
        }
    }
}
void convertqpoly(Polyq *qpoly , Polyp *ppoly){
    for(int j = 0 ; j < 128; j++){
        qpoly->polyarray[j] = ppoly->polyarray[j];
    }
}

void convertqvec(polyvecq *qvec  , polyvecp *pvec){
    for(int i = 0 ; i < qvec->k; i++){
        convertqpoly(qvec->polyarray[i] , pvec->polyarray[i]);
    }
}

void generatem(Polyq **m1list , unsigned char* rseed , polyvecq *x , polyvecq*e2 , unsigned char*H1input){
    
    
    //x || a1 || e2 || a2
    for(int i = 0;i<N;i++){
        m1list[i] = new Polyq(0);
        m1list[i]->copy(x->polyarray[i]);
        m1list[i+1+N] = new Polyq(0);
        m1list[N+1+i]->copy(e2->polyarray[i]);
    }
    
    m1list[N] = new Polyq(0);
    geta(m1list[N] , BETAE - x->norm());
    m1list[N*2+1] = new Polyq(0);
    geta(m1list[N*2+1] ,BETAE - e2->norm());
    
    //tag 
    m1list[N*2+2] = new Polyq(0);
    m1list[N*2+3] = new Polyq(0);
    for(int i = 0 ; i < 16;i++){
        for(int j = 0 ; j < 8 ; j++){
            m1list[N*2+2]->polyarray[i*8+j] = (H1input[i]>>(7-j))&1;
            m1list[N*2+3]->polyarray[i*8+j] = (H1input[i+16]>>(7-j))&1;
        }
    }

    // r || a3

    keccak_state rstate;
    shake256_init(&rstate);
    shake256_absorb(&rstate , rseed , 32);
    int64_t normR = 0;
    for(int i = 0 ; i < K1;i++){
        m1list[N*2+4+i] = new Polyq(&rstate , 3 , 0);
        normR += m1list[N*2+4+i]->norm();
    }
    m1list[N*2+4+K1] = new Polyq(0);
    geta(m1list[N*2+4+K1] ,BETAR - normR);
}
void computeV(polyvecq *v0 ,polyvecq *v1,unsigned char*pke ,unsigned char* rseed ,polyvecq *x ,unsigned char* H1input ,unsigned char* uv){
    polymatp *A = new polymatp(pke,32,0,4,4,1);
    polymatq *C = new polymatq(4,4);

    //convert q
    for(int i = 0;i<4;i++){
        for(int j = 0 ; j < 4;j++){
            for(int k = 0 ; k < 128;k++){
                C->vecarray[i]->polyarray[j]->polyarray[k] = A->vecarray[i]->polyarray[j]->polyarray[k];
            }
            C->vecarray[i]->polyarray[j]->nttflag = 1;
        }
    }

    keccak_state pke_state , rstate;
    shake256_init(&pke_state);
    shake256_absorb(&pke_state , pke+32 , 32);

    shake256_init(&rstate);
    shake256_absorb(&rstate , rseed , 32);
    
    // compute v0
    polyvecq *r = new polyvecq(&rstate , 3 , K1 , 0);
    polyvecq *u = new polyvecq(K1);
    polyvecp *temppoly;
    temppoly = new polyvecp(K1);
    temppoly->from_char(uv , 0 , 0);
    convertqvec(u , temppoly);
    C->right_mul(v0,r);
    v0->to_poly();
    v0->sub(v0 , u);
    int64_t temp;
    for(int i = 0; i < K1;i++){
        for(int j = 0;j<128;j++){
            temp = v0->polyarray[i]->polyarray[j];
            temp = (int)((double)temp/3329 + 0.5);
            v0->polyarray[i]->polyarray[j] = temp;
        }
    }

    //compute v1
    polyvecq *t;
    
    Polyp *v;
    unsigned char mchar[32 + VECETALEN];
    for(int i = 0 ; i < K3 ;i++){
        temppoly = new polyvecp(&pke_state , 0 , K1 , 1);
        t = new polyvecq(K1);
        convertqvec(t , temppoly);
        t->mul(v1->polyarray[i] , r);
        encodem(v1->polyarray[i] , mchar+i*16);
        v = new Polyp(0);
        v->from_char(uv+POLYPLEN*(i+1),0,0);
        for(int j = 0 ; j < 128 ;j++){
            v1->polyarray[i]->polyarray[j] = (int)((double)(v1->polyarray[i]->polyarray[j] - v->polyarray[j])/3329+0.5);
        }
    }
}
void computeh(polyvecq* hvec ,polyvecq*g,Polyq *F,Polyq *G[],polymatq *Gamma){
    Polyq temp(1);
    hvec->to_ntt();
    for(int j = 0 ; j < TAU ; j++){
        Gamma->vecarray[j]->polyarray[N]->mul(&temp , F);
        temp.add(hvec->polyarray[j] , hvec->polyarray[j]);
        for(int i = 0 ; i < 3 ; i++){
            Gamma->vecarray[j]->polyarray[N+i+1]->mul(&temp , G[i]);
            temp.add(hvec->polyarray[j] , hvec->polyarray[j]);
        }
    }
    hvec->to_poly();
    hvec->add(hvec , g);
}
void D2mul(Polyq *res , Polyq* D2list[] , polyvecq *left , polyvecq *right){
    Polyq temp(1);
    int Len = 2*N+M+K1+3;
    res->reset(1);
    for(int i = 0 ; i < N+1 ; i++){
        left->polyarray[Len+i]->mul(&temp , D2list[0]);
        right->polyarray[i]->mul(&temp , &temp);
        res->add(res , &temp);
    }
    for(int i = N+1 ; i < 2*N+2 ; i++){
        left->polyarray[Len+i]->mul(&temp , D2list[1]);
        right->polyarray[i]->mul(&temp , &temp);
        res->add(res , &temp);
    }
    for(int i = 2*N+2 ; i < 2*N+2+M ; i++){
        left->polyarray[Len+i]->mul(&temp , D2list[2]);

        right->polyarray[i]->mul(&temp , &temp);
        
        res->add(res , &temp);
    }
    for(int i = 2*N+2+M ; i < Len ; i++){
        left->polyarray[Len+i]->mul(&temp , D2list[3]);
        right->polyarray[i]->mul(&temp , &temp);
        res->add(res , &temp);
    }
}
int Prove1(unsigned char* pi , unsigned char* crsseed ,unsigned char*rseed,polyvecq *x , polyvecq* e2, unsigned char* H1input, unsigned char* uv ,\
    unsigned char*pke , polymatq *A , polyvecq *c){

    //unpack crs
    keccak_state crsstate;
    shake256_init(&crsstate);
    randseed(crsseed , 32);
    shake256_absorb(&crsstate , crsseed,32);

    polymatq *E1 = new polymatq(&crsstate, 0, N1 ,2*N+M+K1+3 , 1);
    polymatq *E2 = new polymatq(&crsstate, 0, N1 ,M2 , 1);
    polymatq *F1 = new polymatq(&crsstate, 0, N ,M2 , 1);
    polymatq *Fg = new polymatq(&crsstate, 0, TAU ,M2, 1);
    polymatq *F3 = new polymatq(&crsstate, 0, 2 ,M2 , 1);
    polymatq *F4 = new polymatq(&crsstate, 0, 2 ,M2 , 1);
    polymatq *F5 = new polymatq(&crsstate, 0, 3 ,M2 , 1);
    polyvecq *fg = new polyvecq(&crsstate, 0, M2, 1);
    // generate m
    Polyq *m1list[2*N+M+K1+3];
    generatem(m1list , rseed , x , e2 , H1input);
    
    polyvecq *m1 = new polyvecq(m1list , 2*N+M+K1+3);

    keccak_state m2state;
    shake256_init(&m2state);
    shake256_absorb(&m2state , H1input,64);
    
    polyvecq *m2 = new polyvecq(&m2state , 0 , N , 0);
    //generate e
    
reboot:
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
    polyvecq *tF = new polyvecq(N);
    
    E1->right_mul(&temp1 ,m1);
    E2->right_mul(&temp2 , e);
    
    temp1.add(tE , &temp2);
    tE->to_poly();
    polyvecq temp3(N);
    F1->right_mul(&temp3 , e);
    temp3.to_poly();
    temp3.add(tF , m2);

    //FR step 2
    polyvecq *y1 = new polyvecq(Y1 , 2*N+M+K1+3,0);
    polyvecq *y2 = new polyvecq(Y2 , M2,0);
    polyvecq *w = new polyvecq(N1);

    E1->right_mul(&temp1 , y1);
    E2->right_mul(&temp2 , y2);
    temp1.add(w , &temp2);
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
    
    polyvecq *y4 = new polyvecq(Y4 , 2 , 0);
    
    polyvecq *t3 = new polyvecq(2);
    polyvecq *t4 = new polyvecq(2);
    
    F3->right_mul(t3 , e);
    t3->to_poly();
    t3->add(t3 , y3);
    F4->right_mul(t4 , e);
    t4->to_poly();
    t4->add(t4 , y4);

    //FR step 6
    int b[3];
    
    polyvecq *t5 = new polyvecq(3);
    F5->right_mul(t5 , e);
    for(int i = 0 ; i < 3;i++){
        b[i] = 2*(rand()%2)-1;
        t5->polyarray[i]->polyarray[0] =t5->polyarray[0]->mod(t5->polyarray[i]->polyarray[0]+b[i]);
    }

    //SR
    unsigned char Rseed[POLYLEN*(7+TAU+N1*2+N)+32+VECETALEN];
    memcpy(Rseed , crsseed , 32);


    // todo: insert statement x
    memset(Rseed+32 , 0 , VECETALEN);

    t3->to_char(Rseed+32+VECETALEN , 0);
    t4->to_char(Rseed+32+VECETALEN+POLYLEN*2 , 0);
    t5->to_char(Rseed+32+VECETALEN+POLYLEN*4 , 0);
    tg->to_char(Rseed+32+VECETALEN+POLYLEN*7 , 0);
    tE->to_char(Rseed+32+VECETALEN+POLYLEN*(7+TAU) , 0);
    tF->to_char(Rseed+32+VECETALEN+POLYLEN*(7+TAU+N1) , 0);
    w->to_char(Rseed+32+VECETALEN+POLYLEN*(7+TAU+N1+N) , 0);

    keccak_state Rstate;
    shake256_init(&Rstate);
    shake256_absorb(&Rstate , Rseed,POLYLEN*(7+TAU+N1*2+N)+32+VECETALEN);
    
    int64_t *R0 = (int64_t*)malloc(sizeof(int64_t)*256*(2*N+M+K1+3)*128);
    int64_t *R1 = (int64_t*)malloc(sizeof(int64_t)*256*(K2+K3)*128);

    cbd(&Rstate , R0 , 1,256*128*(2*N+M+K1+3));
    cbd(&Rstate , R1 ,1, 256*128*(K2+K3));

    //SR step3
    polyvecq *v0 = new polyvecq(K2);
    polyvecq *v1 = new polyvecq(K3);
    computeV(v0 , v1 , pke , rseed , x , H1input , uv);

    //SR step4
    int64_t *m1vec = (int64_t*)malloc(sizeof(int64_t)*(2*N+M+K1+3)*128);
    int64_t z3vec[256];
    int64_t r0m1[256];
    int64_t y3vec[256];
    polyvec2vec(m1vec , m1);
    polyvec2vec(y3vec , y3);
    matrix_mul(r0m1 , R0 , m1vec , 256 , (2*N+M+K1+3)*128);
    memcpy(z3vec , r0m1 , 256*sizeof(int64_t));
    vec_mul_si(z3vec , b[0] , 256);
    vec_add(z3vec , y3vec , 256);

    int64_t z4vec[256];
    int64_t r1v_[256];
    int64_t v_vec[(K2+K3)*128];
    int64_t y4vec[256];
    
    polyvec2vec(v_vec , v0);
    polyvec2vec(v_vec+K2 , v1);

    polyvec2vec(y4vec , y4);
    matrix_mul(r1v_ , R1 , v_vec , 256 , (K2+K3)*128);
    
    memcpy(z4vec , r1v_ , 256*sizeof(int64_t));
    vec_mul_si(z4vec , b[1] , 256);
    vec_add(z4vec , y4vec , 256);
    
    if(Rej0(z3vec ,r0m1 , Y3 , 7)){
        return 0;
    }
    if(Rej0(z4vec , r1v_ , Y4 , 1.3)){
        return 0;
    }



    // TR
    unsigned char Gammaseed[256*16];
    vec2char(Gammaseed , z3vec , 256);
    vec2char(Gammaseed+256*8 , z4vec , 256);

    keccak_state gamma_state;
    shake256_init(&gamma_state);
    shake256_absorb(&gamma_state , Rseed,POLYLEN*(7+TAU+N1*2+N)+32+VECETALEN);
    shake256_absorb(&gamma_state , Gammaseed , 256*16);
    
    polymatq *Gamma = new polymatq(TAU , N+4);
    sample_gamma(&gamma_state , Gamma);

    // compute EFG
    //Polyq *Elist[N];
    Polyq *F = new Polyq(0);
    Polyq *G[3];

    //E
    /*
    polyvecq *temp4 = new polyvecq(N);
    polyvecq *temp5 = new polyvecq(N);
    e2->add(temp4 , m2);
    c->to_poly();
    temp4->sub(temp4 , c);
    for(int i = 0 ; i < N ; i++){
        Elist[i] = new Polyq(1);
        A->vecarray[i]->sigma(temp5);
        temp5->mul(Elist[i] , x);
        Elist[i]->to_poly();
        
        Elist[i]->add(Elist[i] , temp4->polyarray[i]);
    }*/
    
    //F 
    Polyq *temp6 = new Polyq(0);
    Polyq *temp7 = new Polyq(0);
    for(int i = 0 ; i < 2;i++){
        temp7->reset(0);
        for(int j = 0 ; j < 128 ; j++){
            temp7->polyarray[j] = 1-m1list[2*N+2+i]->polyarray[j];
        }
        temp7->sigma(temp7);
        temp7->mul(temp6 ,m1list[2*N+2+i]);
        temp6->to_poly();
        temp6->add(F , F);
    }

    //G
    polyvecq *temp8 = new polyvecq(N+1);
    polyvecq *temp9 = new polyvecq(N+1);
    for(int i = 0 ; i < N+1 ; i++){
        temp8->polyarray[i]->copy(m1list[i]);
    }
    G[0] = new Polyq(1);
    temp8->sigma(temp9);
    temp8->mul(G[0] , temp9);
    G[0]->to_poly();
    G[0]->polyarray[0] -= BETAE;
    delete temp8;
    
    temp8 = new polyvecq(N+1);
    for(int i = 0 ; i < N+1 ; i++){
        temp8->polyarray[i]->copy(m1list[N+1+i]);
    }
    G[1] = new Polyq(1);
    temp8->sigma(temp9);
    temp8->mul(G[1] , temp9);
    G[1]->to_poly();
    G[1]->polyarray[0] -= BETAE;
    delete temp8;
    delete temp9;

    temp8 = new polyvecq(K1+1);
    temp9 = new polyvecq(K1+1);
    for(int i = 0 ; i < K1+1 ; i++){
        temp8->polyarray[i]->copy(m1list[N*2+M+2+i]);
    }

    G[2] = new Polyq(1);
    temp8->sigma(temp9);
    temp8->mul(G[2] , temp9);
    G[2]->to_poly();
    G[2]->polyarray[0] -= BETAR;
    delete temp8;
    delete temp9;

    //compute hvec
    polyvecq *hvec = new polyvecq(TAU);
    computeh(hvec , g ,F,G,Gamma);

    //FourR step1

    //H3input

    unsigned char museed[TAU*POLYLEN];
    hvec->to_char(museed);
    keccak_state mu_state;
    shake256_init(&mu_state);
    shake256_absorb(&mu_state , Rseed,POLYLEN*(7+TAU+N1*2+N)+32+VECETALEN);
    shake256_absorb(&mu_state , Gammaseed , 256*16);
    shake256_absorb(&mu_state , museed , TAU*POLYLEN);

    polyvecq *mu = new polyvecq(&mu_state ,0, TAU , 0);



    //FourR step 2
    polymatq *Fymat = new polymatq(4+TAU+N , M2);
    polyvecq *y = new polyvecq((2*N+M+K1+3)*2 + (4+TAU+N)*2);
    constructFy(Fymat , F1 , Fg , F3 , F4);
    constructy(y , Fymat , y1 , y2);

    //FourR step 3
    Polyq *D2list[4];
    constructD2(D2list , mu , Gamma);


    polyvecq *d1vec = new polyvecq(2*(2*N+M+K1+3)+2*(N+TAU+4));
    constructd1(d1vec , mu , Gamma , A);

    //FourR step 4

    //construct m
    polyvecq *m = new polyvecq(2*(2*N+M+K1+3)+2*(N+TAU+4));
    constructm(m , m1 , m2 , g , y3 , y4);




    //compute g1
    Polyq tempg(1);
    Polyq *g1 = new Polyq(1);
    D2mul(&tempg , D2list , m , y);
    g1->add(g1 , &tempg);
    D2mul(&tempg , D2list,y , m);
    g1->add(g1 , &tempg);
    d1vec->mul(&tempg , y);
    g1->add(g1 , &tempg);

    //compute t
    Polyq *t = new Polyq(1);
    fg->mul(t , e);
    
    t->add(t , g1);


    //compute v
    Polyq *v = new Polyq(1);
    D2mul(&tempg , D2list , y ,y);
    v->add(v , &tempg);
    fg->mul(&tempg , y2);
    v->add(v , &tempg);


    //FifthR

    //H4input
    unsigned char cseed[POLYLEN*2];

    t->to_char(cseed);
    v->to_char(cseed+POLYLEN);
    keccak_state c_state;
    shake256_init(&c_state);
    shake256_absorb(&c_state , Rseed,POLYLEN*(7+TAU+N1*2+N)+32+VECETALEN);
    shake256_absorb(&c_state , Gammaseed , 256*16);
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
    //cpoly->mul_num(cpoly , b[2]);
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
    
    if(z1->norm() > B1*B1 || z2->norm() > B2*B2 || vec_norm(z3vec , 256) > B3*B3 || vec_norm(z4vec , 256) > B4*B4){
        return 0;
    }
    memcpy(pi , Rseed+VECETALEN+32 , POLYLEN*(7+TAU+N1*2+N));//POLYLEN*(7+TAU+N1*2+N)
    memcpy(pi+POLYLEN*(7+TAU+N1*2+N) , Gammaseed,256*16);
    memcpy(pi + POLYLEN*(7+TAU+N1*2+N) + 256*16 , museed, TAU*POLYLEN);
    memcpy(pi + POLYLEN*(7+2*TAU+N1*2+N) + 256*16 , cseed , POLYLEN*2);
//2*N+M+K1+3
    z1->to_char(pi + POLYLEN*(9+2*TAU+N1*2+N) + 256*16);
    z2->to_char(pi + POLYLEN*(12+2*TAU+N1*2+N*3+M+K1)+256*16);
    // 12+2*TAU+N1*2+N*3+M+K1+M2

    
    return 1;
    
}

int Verify1(unsigned char *crsseed, unsigned char *pi , polymatq *A, polyvecq *c){
    keccak_state crsstate;
    shake256_init(&crsstate);
    shake256_absorb(&crsstate , crsseed,32);

    polymatq *E1 = new polymatq(&crsstate, 0, N1 ,2*N+M+K1+3 , 1);
    polymatq *E2 = new polymatq(&crsstate, 0, N1 ,M2 , 1);
    polymatq *F1 = new polymatq(&crsstate, 0, N ,M2 , 1);
    polymatq *Fg = new polymatq(&crsstate, 0, TAU ,M2, 1);
    polymatq *F3 = new polymatq(&crsstate, 0, 2 ,M2 , 1);
    polymatq *F4 = new polymatq(&crsstate, 0, 2 ,M2 , 1);
    polymatq *F5 = new polymatq(&crsstate, 0, 3 ,M2 , 1);
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
    polyvecq *tF = new polyvecq(N);
    polyvecq *w = new polyvecq(N1);
    t3->from_char(Rseed+32+VECETALEN , 0);
    t4->from_char(Rseed+32+VECETALEN+POLYLEN*2 , 0);
    tg->from_char(Rseed+32+VECETALEN+POLYLEN*7 , 0);
    tE->from_char(Rseed+32+VECETALEN+POLYLEN*(7+TAU) , 0);
    tF->from_char(Rseed+32+VECETALEN+POLYLEN*(7+TAU+N1) , 0);
    w->from_char(Rseed+32+VECETALEN+POLYLEN*(7+TAU+N1+N) , 1);
    keccak_state Rstate;
    shake256_init(&Rstate);
    shake256_absorb(&Rstate , Rseed,POLYLEN*(7+TAU+N1*2+N)+32+VECETALEN);
    
    int64_t *R0 = (int64_t*)malloc(sizeof(int64_t)*256*(2*N+M+K1+3)*128);
    int64_t *R1 = (int64_t*)malloc(sizeof(int64_t)*256*(K1+K2)*128);

    cbd(&Rstate , R0 , 1,256*128*(2*N+M+K1+3));
    cbd(&Rstate , R1 , 1 , 256*128*(K2+K3));


    //unpack alpha2

    int64_t z3vec[256];
    int64_t z4vec[256];
    unsigned char Gammaseed[256*16];
    memcpy(Gammaseed,pi+POLYLEN*(7+TAU+N1*2+N) ,256*16);
    char2vec(z3vec , Gammaseed , 256);
    char2vec(z4vec , Gammaseed+256*8 , 256);

    keccak_state gamma_state;
    shake256_init(&gamma_state);
    shake256_absorb(&gamma_state , Rseed,POLYLEN*(7+TAU+N1*2+N)+32+VECETALEN);
    shake256_absorb(&gamma_state , Gammaseed , 256*16);
    polymatq *Gamma = new polymatq(TAU , N+4);
    sample_gamma(&gamma_state , Gamma);

    // unpack alpha3
    polyvecq *hvec = new polyvecq(TAU);

    unsigned char museed[TAU*POLYLEN];
    memcpy(museed, pi + POLYLEN*(7+TAU+N1*2+N) + 256*16 ,TAU*POLYLEN);
    hvec->from_char(museed , 0);
    keccak_state mu_state;
    shake256_init(&mu_state);
    shake256_absorb(&mu_state , Rseed,POLYLEN*(7+TAU+N1*2+N)+32+VECETALEN);
    shake256_absorb(&mu_state , Gammaseed , 256*16);
    shake256_absorb(&mu_state , museed , TAU*POLYLEN);

    
    polyvecq *mu = new polyvecq(&mu_state ,0, TAU , 0);

    //unpack alpha4
    unsigned char cseed[POLYLEN*2];
    memcpy(cseed ,pi + POLYLEN*(7+2*TAU+N1*2+N) + 256*16 , POLYLEN*2);
    Polyq *t = new Polyq(1);
    Polyq *v = new Polyq(1);
    t->from_char(cseed , 1);
    v->from_char(cseed+POLYLEN , 1);
    keccak_state c_state;
    shake256_init(&c_state);
    shake256_absorb(&c_state , Rseed,POLYLEN*(7+TAU+N1*2+N)+32+VECETALEN);
    shake256_absorb(&c_state , Gammaseed , 256*16);
    shake256_absorb(&c_state , museed , TAU*POLYLEN);
    shake256_absorb(&c_state , cseed , POLYLEN*2);
    Polyq *cpoly = new Polyq(0);
    Csamplepoly(&c_state , cpoly);

    //constructz
    polyvecq *z1 = new polyvecq(2*N+M+K1+3);
    polyvecq *z2 = new polyvecq(M2);
    polyvecq *z = new polyvecq(2*(2*N+M+K1+3)+2*(N+TAU+4));
    z1->from_char(pi + POLYLEN*(9+2*TAU+N1*2+N) + 256*16 , 0);
    z2->from_char(pi + POLYLEN*(12+2*TAU+N1*2+N*3+M+K1)+256*16 , 0);

    polymatq *Fymat = new polymatq(4+TAU+N , M2);
    constructFy(Fymat , F1 , Fg , F3 , F4);
    constructz(z , Fymat , tF , tg , t3,t4,cpoly , z1 , z2);



    // constructD
    Polyq *D2list[4];
    constructD2(D2list , mu , Gamma);

    polyvecq *d1vec = new polyvecq(2*(2*N+M+K1+3)+2*(N+TAU+4));
    constructd1(d1vec , mu , Gamma , A);

    Polyq *d0 = new Polyq(0);
    constructd0(d0 , mu , Gamma , c , hvec);

    //check norm
    
    if(z1->norm() > B1*B1 || z2->norm() > B2*B2 || vec_norm(z3vec , 256) > B3*B3 || vec_norm(z4vec , 256) > B4*B4){
        std::cout << "wrong0" << std::endl;
        return 0;
    }
    
    // check E1z1 + E2z2 == w+ ctE

    polyvecq left(N1 , 1) , temp(N1 , 1) , right(N1 , 1);
    E1->right_mul(&temp , z1);
    left.add(&left, &temp);
    E2->right_mul(&temp , z2);
    left.add(&left , &temp);
    left.to_poly();

    tE->mul_poly(&right , cpoly);
    right.add(&right , w);
    
    right.to_poly();
    if(left.equal(&right) == 0){
        std::cout << "wrong1" << std::endl;
        return 0;
    }

    // check z^TD2z + cd1^Tz + c^2d0 + fg^Tz2= v + ct

    Polyq leftpoly(1) ,temppoly(1);
    
    D2mul(&leftpoly , D2list ,z ,  z);
    d1vec->mul(&temppoly , z);
    temppoly.mul(&temppoly , cpoly);
    leftpoly.add(&leftpoly , &temppoly);

    d0->mul(&temppoly , cpoly);
    cpoly->mul(&temppoly , &temppoly);
    leftpoly.add(&leftpoly , &temppoly);
    
    fg->mul(&temppoly , z2);
    leftpoly.add(&leftpoly , &temppoly);

    cpoly->mul(&temppoly , t);
    leftpoly.sub(&leftpoly , &temppoly);

    if(leftpoly.equal(v) == 0){
        std::cout << "wrong3" << std::endl;
        return 0;
    }

    hvec->to_poly();
    //check h
    for(int i = 0 ; i < TAU ;i++){
        if(hvec->polyarray[i]->polyarray[0] !=0 ){
            std::cout << "wrong4" << std::endl;
            return 0;
        }
    }
    return 1;
}
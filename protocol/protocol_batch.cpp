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

void generate_batch_m1(Polyq *m1list[] ,int b[], polyvecq* s0, polyvecq* s1, polyvecq* e0,\
polyvecq* e1, polyvecq *e3){
//m1list 4N+N_BATCH+6
    s0->to_poly();
    s1->to_poly();
    e0->to_poly();
    e1->to_poly();
    m1list[0] = new Polyq(0);
    for(int i=0; i < N_BATCH; i++){
        m1list[0]->polyarray[i] = b[i];
    }
    
// s0 a0 s1 a1 e0 a2 e1 a3
    for(int i = 0 ; i<N ;i++){
        m1list[i+1] = new Polyq(0);
        m1list[i+1]->copy(s0->polyarray[i]);
        m1list[i+1+N+1] = new Polyq(0);
        m1list[i+1+N+1]->copy(s1->polyarray[i]);
        m1list[i+1+2*(N+1)] = new Polyq(0);
        m1list[i+1+2*(N+1)]->copy(e0->polyarray[i]);
        m1list[i+1+3*(N+1)] = new Polyq(0);
        m1list[i+1+3*(N+1)]->copy(e1->polyarray[i]);
    }

//a0 a1 a2 a3    
    m1list[N+1] = new Polyq(0);
    geta(m1list[N+1] , BETAE - s0->norm());
    m1list[N*2+2] = new Polyq(0);
    geta(m1list[N*2+2] ,BETAE - s1->norm());
    m1list[N*3+3] = new Polyq(0);
    geta(m1list[N*3+3] , BETAE - e0->norm());
    m1list[N*4+4] = new Polyq(0);
    geta(m1list[N*4+4] ,BETAE - e1->norm());
//e_3
    for(int i = 0;i<N_BATCH;i++){
        m1list[i+N*4+5] = new Polyq(0);
        m1list[i+N*4+5]->copy(e3->polyarray[i]);
        // std::cout<<"e3 "<<i<<"-th entry:"<<std::endl;
        // e3->polyarray[i]->output();
    }
//a_4
    m1list[N*4+5+N_BATCH] = new Polyq(0);
    geta(m1list[N*4+5+N_BATCH] , NBETAE2-e3->norm());
}


void constructFy_batch(polymatq *Fymat ,polymatq *Fg ,polymatq *F3){

    for(int i = 0 ; i < TAU ; i++){
        Fymat->vecarray[i] = Fg->vecarray[i];
    }
    for(int i = 0 ; i < 2;i++){
        Fymat->vecarray[TAU+i] = F3->vecarray[i];
    }
}
void constructy_batch(polyvecq *y ,polymatq *Fymat ,polyvecq *y1 ,polyvecq*y2){
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
void computeh(polyvecq* hvec, polyvecq*g, Polyq *F[],polymatq *Gamma){
    Polyq temp(1);
    hvec->to_ntt();
    for(int j = 0 ; j < TAU ; j++){
        for(int i = 0 ; i < 5 ; i++){
            Gamma->vecarray[j]->polyarray[2*N+N_BATCH+i]->mul(&temp , F[i]);
            temp.add(hvec->polyarray[j] , hvec->polyarray[j]);
        }
    }
    hvec->to_poly();
    //hvec->add(hvec , g);
}

//D2 (coefficients)
void constructD2_batch(Polyq *D2list[] , polyvecq *mu , polymatq *Gamma, polymatq *c){

    polyvecq temp1(4);
        for(int j = 0 ; j < TAU ; j++){
            mu->polyarray[j]->mul(temp1.polyarray[0] , Gamma->vecarray[j]->polyarray[N_BATCH+2*N+5]);
            mu->polyarray[j]->mul(temp1.polyarray[1] , Gamma->vecarray[j]->polyarray[N_BATCH+2*N]);
            mu->polyarray[j]->mul(temp1.polyarray[2] , Gamma->vecarray[j]->polyarray[N_BATCH+2*N+1]);
            mu->polyarray[j]->mul(temp1.polyarray[3] , Gamma->vecarray[j]->polyarray[N+5]);  
        }
        for(int i=0;i<4;i++){
            D2list[i]->copy(temp1.polyarray[i]);//ntt
        }

        // mu_j gamma c_i
        polyvecq temp2(N,1), temp3(N,1), temp4(N,1),temp5(N,1);
        for(int j = 0 ; j < TAU ; j++){
            temp2.reset(1);
            for (int i = 0; i < N_BATCH; i++){
                c->vecarray[i]->mul_poly(&temp5, Gamma->vecarray[j]->polyarray[i]);
                temp2.add(&temp2,&temp5);
            }
            temp2.mul_poly(&temp3, mu->polyarray[j]);
            temp4.add(&temp4, &temp3);
        }
        for(int i=0; i<N;i++){
            D2list[i+4]->copy(temp4.polyarray[i]);//ntt
        }
        
}

//d1
void constructd1_batch(polyvecq *d1vec , polyvecq* mu , polymatq *Gamma ,polymatq *A , polymatq *c){
    
    Polyq temp1(1) , temp2(1) ;
    polyvecq temp3(N), temp8(N,1), temp9(N,1), temp10(N , 1);
    polyvecq temp4(N) , temp5(N , 1), temp6(N_BATCH),temp7(N_BATCH);

     d1vec->polyarray[0]->reset(1);
// int vec b
    for(int j = 0 ; j < TAU ; j++){
        mu->polyarray[j]->mul(&temp1, Gamma->vecarray[j]->polyarray[N_BATCH+2*N+5]);
        //temp1 ntt, temp2 ntt
        temp2.add(&temp2, &temp1);
    }
    temp2.to_ntt();
    d1vec->polyarray[0]->copy(&temp2);
    temp1.reset(1);
    temp2.reset(1);

//s0 
//s1
    A->trans(); //transform
    for(int j = 0 ; j < TAU ; j++){
        temp4.reset(1);
        for(int i = 0 ; i < N ; i++){
            A->vecarray[i]->mul_poly(&temp3 , Gamma->vecarray[j]->polyarray[i+N_BATCH]);
            A->vecarray[i]->mul_poly(&temp8 , Gamma->vecarray[j]->polyarray[i+N_BATCH+N]);
            temp4.add(&temp4 , &temp3);//all ntt
            temp9.add(&temp9 , &temp8); //all ntt
        }
        for(int i = 0 ; i < N_BATCH ; i++){
            c->vecarray[i]->mul_poly(&temp3 , Gamma->vecarray[j]->polyarray[i]);
            temp4.add(&temp4 , &temp3);
            temp9.add(&temp9 , &temp3);
        }      
        temp4.mul_poly(&temp4 , mu->polyarray[j]);
        temp9.mul_poly(&temp9 , mu->polyarray[j]);
        temp5.add(&temp5 , &temp4);
        temp10.add(&temp10 , &temp9);
    }
    A->trans(); //transform back
    temp5.to_ntt();
    temp10.to_ntt();
    for(int i = 0 ; i < N ; i++){
        d1vec->polyarray[i+1]->copy(temp5.polyarray[i]);
        d1vec->polyarray[i+N+2]->copy(temp5.polyarray[i]);
    }

//e0
//e1
    temp5.reset(1);
    temp10.reset(1);
    for(int j = 0 ; j < TAU ; j++){
        for(int i = 0 ; i < N;i++){
            temp4.polyarray[i]->copy(Gamma->vecarray[j]->polyarray[i+N_BATCH]);
            temp9.polyarray[i]->copy(Gamma->vecarray[j]->polyarray[i+N_BATCH+N]);
        }
        temp4.mul_poly(&temp4 , mu->polyarray[j]);
        temp9.mul_poly(&temp9 , mu->polyarray[j]);
        temp5.add(&temp5 , &temp4);
        temp10.add(&temp10 , &temp9);
    }
    temp5.to_ntt();
    temp10.to_ntt();
    for(int i = 0 ; i < N ; i++){
        d1vec->polyarray[i+2*N+3]->copy(temp5.polyarray[i]);
        d1vec->polyarray[i+3*N+4]->copy(temp10.polyarray[i]);
    }
    temp5.reset(1);
    temp10.reset(1);

//vec e3 length N_BATCH
    temp7.reset(1);
    for(int j = 0 ; j < TAU ; j++){ 
        for(int i = 0 ; i < N_BATCH; i++){
            temp6.polyarray[i]->copy(Gamma->vecarray[j]->polyarray[i]);
        }
        temp6.mul_poly(&temp6 , mu->polyarray[j]);
        temp7.add(&temp7 , &temp6);
    }
    temp7.to_ntt();
    for(int i = 0 ; i < N_BATCH ; i++){
        d1vec->polyarray[i+4*N+5]->copy(temp7.polyarray[i]);
    }
}

//m = m_1 || \sigma(m_1) || g ||y_3 || \sigma(g||y_3) 
void constructm_batch(polyvecq *m ,polyvecq *m1 ,polyvecq *g ,polyvecq *y3){
    int mlen = M1_BATCH;
    int m2len = TAU+2;
    for(int i = 0 ; i < mlen; i++){
        m->polyarray[i]->copy(m1->polyarray[i]);
        m1->polyarray[i]->sigma(m->polyarray[i+mlen]);
        // m->polyarray[i]->output();
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

//z
void constructz_batch(polyvecq *z ,polymatq*Fymat, polyvecq *tg, polyvecq *t3,Polyq*cpoly, polyvecq *z1, polyvecq *z2){
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

// result of left D_2 right
void D2mul_batch(Polyq *res , Polyq* D2list[] , polyvecq *left , polyvecq *right){
    Polyq temp(1), temp1(1), temp2(1),temp3(1), temp4(1);
    int Len = M1_BATCH;
    int Len1 = N+1;
    res->reset(1); //ntt
    
    //left vec b
    left->polyarray[Len]->mul(&temp , D2list[0]);
    right->polyarray[0]->mul(&temp , &temp);
    res->sub(res , &temp); //ntt

    //left * c1,...cN *vec b
    temp3.reset(1);
    temp4.reset(1);
    for(int i= 0; i<N ; i++){
        left->polyarray[Len+1+i]->mul(&temp1, D2list[i+4]); //mul(s0, ci)
        left->polyarray[Len+1+i+N+1]->mul(&temp2, D2list[i+4]); //mul(s1, ci)
        temp1.mul_num(&temp1, right->polyarray[0]->polyarray[i]); //mul(s0, ci, bi)
        temp2.mul_num(&temp2, right->polyarray[0]->polyarray[i]); //mul(s1, ci, bi)
        temp1.to_ntt();
        temp2.to_ntt();
        temp3.add(&temp3, &temp1); //add mul(s0, ci, bi)  //ntt
        temp4.add(&temp4, &temp2); //add mul(s0, ci, bi)  //ntt
    }

    for(int i = 1 ; i < 1+Len1 ; i++){
        left->polyarray[Len+i]->mul(&temp , D2list[1]);
        right->polyarray[i]->mul(&temp , &temp);
        res->add(res , &temp);
    }
    for(int i = 1+Len1 ; i < 1+2*Len1 ; i++){
        left->polyarray[Len+i]->mul(&temp , D2list[2]);
        right->polyarray[i]->mul(&temp , &temp);
        res->add(res , &temp);
    }
    for(int i = 2+2*Len1 ; i < 2 ; i++){
        left->polyarray[Len+i]->mul(&temp , D2list[3]);
        right->polyarray[i]->mul(&temp , &temp);
        res->add(res , &temp);
    }
}

//d0
void constructd0_batch(Polyq *d, polyvecq *mu, polymatq *Gamma, polyvecq*h , polyvecq* hvec, polyvecq*t0, polyvecq*t1){
    Polyq temp1(1) , temp2(1);
    Polyq temp3(1) , temp4(1);
    d->to_ntt();
    mu->to_ntt();
    hvec->mul(d , mu);//mu_j h, d:ntt

//mu_j gamma hi
    temp2.reset(1);
    temp1.reset(1);
    temp3.reset(1);
    for(int j = 0 ; j < TAU ; j++){
        for(int i=0;i<N_BATCH;i++){
            h->polyarray[i]-> mul(&temp1 , Gamma->vecarray[j]->polyarray[i]); //gamma_{j,i} hi
            temp1.to_ntt();
            temp2.add(&temp2 , &temp1);//ntt
        }
        mu->polyarray[j]->mul(&temp3 , &temp2); //mu_j sum gamma hi
    }
    temp3.to_ntt();    
    d->add(d , &temp3); //ntt form

    //mu_j gamma t0
    temp2.reset(1);
    temp1.reset(1);
    temp3.reset(1);
    for(int j = 0 ; j < TAU ; j++){
        for(int i=0;i<N;i++){
            t0->polyarray[i]-> mul(&temp1 , Gamma->vecarray[j]->polyarray[i+N_BATCH]); //gamma_{j,i} t0
            temp1.to_ntt();
            temp2.to_ntt();
            temp2.add(&temp2 , &temp1);
        }
        mu->polyarray[j]->mul(&temp3 , &temp2); //mu_j sum gamma t0
    }
    d->add(d , &temp3); //ntt form
        
//mu_j gamma t1
    temp2.reset(1);
    temp1.reset(1);
    temp3.reset(1);
    for(int j = 0 ; j < TAU ; j++){
        for(int i=0;i<N;i++){
            t1->polyarray[i]-> mul(&temp1 , Gamma->vecarray[j]->polyarray[i+N_BATCH+N]); //gamma_{j,i} t1
            temp1.to_ntt();
            temp2.to_ntt();
            temp2.add(&temp2 , &temp1);
        }
        mu->polyarray[j]->mul(&temp3 , &temp2); //mu_j sum gamma t1
    }
    d->add(d , &temp3); //ntt
 
 
 //beta_r   
    temp2.reset(0);
    temp1.reset(0);
    temp3.reset(0);
    for(int j = 0 ; j < TAU ; j++){
        for(int i=0;i<4;i++){
            temp1.polyarray[0] = BETAE * Gamma->vecarray[j]->polyarray[i+N_BATCH+2*N]->polyarray[0]; //gamma beta_r 
            temp1.to_ntt();
            temp2.to_ntt();
            temp2.add(&temp2 , &temp1);
        }
        mu->polyarray[j]->mul(&temp3 , &temp2); //mu_j sum gamma beta_r
    }
    d->add(d , &temp3); //ntt

    //beta_r 
    temp2.reset(0);
    temp1.reset(0);
    temp3.reset(0);
    for(int j = 0 ; j < TAU ; j++){
        temp1.polyarray[0] = (N_BATCH) * BETAE2 * Gamma->vecarray[j]->polyarray[N_BATCH+2*N+4]->polyarray[0]; //gamma beta_e
        mu->polyarray[j]->mul(&temp2 , &temp1);
    }
    d->add(d , &temp2); //ntt
    d->to_poly();
    d->mul_num(d , -1);
}


int Prove_batch(unsigned char* pi , unsigned char* crsseed , int b[], polyvecq *s0, polyvecq *s1, polyvecq *e0,  polyvecq *e1,\
    polyvecq *e3, polyvecq* h, polymatq*A, polyvecq*t0, polyvecq *t1, polymatq *c){
    //unpack crs
    keccak_state crsstate;
    shake256_init(&crsstate);
    randseed(crsseed , 32);
    shake256_absorb(&crsstate , crsseed, 32);

    polymatq *E1 = new polymatq(&crsstate, 0, N1 , M1_BATCH , 1);
    polymatq *E2 = new polymatq(&crsstate, 0, N1 ,M2 , 1);
    polymatq *Fg = new polymatq(&crsstate, 0, TAU ,M2, 1);
    polymatq *F3 = new polymatq(&crsstate, 0, 2 ,M2 , 1);
    polymatq *F4 = new polymatq(&crsstate, 0, 2 ,M2 , 1);
    polyvecq *fg = new polyvecq(&crsstate, 0, M2, 1);
    // generate m

    // std::cout << "E1 in PV: " << std::endl;
    // E1->vecarray[0]->polyarray[0]->output();
    // std::cout << "E2 in PV: " << std::endl;
    // E2->vecarray[0]->polyarray[0]->output();
    
    Polyq *m1list[M1_BATCH];
    generate_batch_m1(m1list , b, s0, s1, e0, e1, e3);
    
    polyvecq *m1 = new polyvecq(m1list , M1_BATCH); //m1 = m1list

    //generate e
    
    unsigned char random_seed[32];
    keccak_state random_state;
    shake256_init(&random_state);
    randseed(random_seed , 32);
    shake256_absorb(&random_state , random_seed,32);
    
    polyvecq *e = new polyvecq(&random_state , ETAR , M2 , 0); //e len: M2


    // First Round step 1

    polyvecq temp1(N1);
    polyvecq temp2(N1);
    polyvecq *tE = new polyvecq(N1);
    
    E1->right_mul(&temp1 ,m1);
    E2->right_mul(&temp2 , e);
    
    temp1.add(tE , &temp2);
    tE->to_poly();
    // std::cout << "tE in PV: " << std::endl;
    // tE->polyarray[0]->output();

    //First Round  step 2
    polyvecq *y1 = new polyvecq(Y1 , M1_BATCH, 0);
    polyvecq *y2 = new polyvecq(Y2 , M2, 0);
    polyvecq *w = new polyvecq(N1);

    E1->right_mul(&temp1 , y1);
    E2->right_mul(&temp2 , y2);
    temp1.add(w , &temp2);
    w->to_poly();
    
    // std::cout << "w in PV: " << std::endl;
    // w->polyarray[0]->output();

    //First Round step 3
    polyvecq *g = new polyvecq(&random_state,0 ,TAU , 0);
    polyvecq *tg = new polyvecq(TAU);
    for(int i = 0;i< TAU;i++){
        g->polyarray[i]->polyarray[0] = 0;
    }
    Fg->right_mul(tg , e);
    tg->to_poly();
    tg->add(tg , g);

    //First Round step 4-5
    polyvecq *y3 = new polyvecq(Y3 , 2 , 0);
    
    polyvecq *t3 = new polyvecq(2);
    //polyvecq *t4 = new polyvecq(2);
    
    F3->right_mul(t3 , e);
    t3->to_poly();
    t3->add(t3 , y3);
    
    //First Round step 6
    //F4->right_mul(t4 , e);
   // t4->to_poly();

    // int bits[2];
    // for(int i = 0 ; i < 2;i++){
    //     bits[i] = 2*(rand()%2)-1;
    //     t4->polyarray[i]->polyarray[0] =t4->polyarray[0]->mod(t4->polyarray[i]->polyarray[0]+bits[i]);
    // }
    
    //Second R
    unsigned char* Rseed;
    Rseed = (unsigned char*)malloc(POLYLEN*(2+TAU+2*N1)+32+VECETALEN);//VECETALEN = XSTATEMENTLEN
    memcpy(Rseed , crsseed , 32); //crsseed len 32

    // todo: insert statement x
    memset(Rseed+32 , 0, VECETALEN);  //x length VECETALEN

    t3->to_char(Rseed+32+VECETALEN , 0);
    //t4->to_char(Rseed+32+VECETALEN+POLYLEN*2 , 0);
    tg->to_char(Rseed+32+VECETALEN+POLYLEN*2 , 0);
    tE->to_char(Rseed+32+VECETALEN+POLYLEN*(2+TAU) , 0);
    w->to_char(Rseed+32+VECETALEN+POLYLEN*(2+TAU+N1) , 0);

    keccak_state Rstate;
    shake256_init(&Rstate);
    shake256_absorb(&Rstate , Rseed, POLYLEN*(2+TAU+2*N1)+32+VECETALEN); //crs + x + alpha1 
    
    int64_t *R = (int64_t*)malloc(sizeof(int64_t)*256*(M1_BATCH)*128);

    cbd(&Rstate , R , 1 , 256*128*(M1_BATCH));

    //SR step4
    int64_t *m1vec = (int64_t*)malloc(sizeof(int64_t)*(M1_BATCH)*128);
    int64_t *z3vec = (int64_t*)malloc(sizeof(int64_t)*256);  //z3 length 256
    int64_t *Rm1 = (int64_t*)malloc(sizeof(int64_t)*256);
    int64_t *y3vec = (int64_t*)malloc(sizeof(int64_t)*256);
    polyvec2vec(m1vec , m1);
    polyvec2vec(y3vec , y3);
    matrix_mul(Rm1 , R , m1vec , 256 , (M1_BATCH)*128);
    memcpy(z3vec , Rm1 , 256*sizeof(int64_t));
    // vec_mul_si(z3vec , bits[0] , 256);
    vec_add(z3vec , y3vec , 256);
    
    if(Rej0(z3vec ,Rm1 , Y3 , 7)){
        return 0;
    }
    // std::cout << "z3vec in PV: " << std::endl;
    
    // Third Round
    unsigned char Gammaseed[256*8];
    vec2char(Gammaseed , z3vec , 256);

    keccak_state gamma_state;
    shake256_init(&gamma_state);
    shake256_absorb(&gamma_state , Rseed, POLYLEN*(4+TAU+2*N1)+32+VECETALEN);  //crs + x + alpha1 
    shake256_absorb(&gamma_state , Gammaseed , 256*8);  //alpha2
    
    polymatq *Gamma = new polymatq(TAU , N_BATCH+2*N+6);
    sample_gamma(&gamma_state , Gamma);

    // compute H1,...H6
    Polyq *F[5];
    polyvecq *temp8 = new polyvecq(N+1);
    polyvecq *temp9 = new polyvecq(N+1);
    int tempbetalist[2] = {BETAE, NBETAE2};
     
    for(int j = 0 ; j < 4 ; j++){
        F[j] = new Polyq(0);
        for(int i = 0 ; i < N+1 ; i++){
            temp8->polyarray[i]->copy(m1list[1+i+(N+1)*j]);
        }
        temp8->sigma(temp9);
        temp8->mul(F[j] , temp9);
        F[j]->to_poly();
        F[j]->polyarray[0] -= tempbetalist[0];
    }

    F[4] = new Polyq(0);
    polyvecq *temp8_1 = new polyvecq(N_BATCH+1);
    polyvecq *temp9_1 = new polyvecq(N_BATCH+1);
    for(int i = 0 ; i < N_BATCH+1 ; i++){
        temp8_1->polyarray[i]->copy(m1list[4*N+5+i]);
        // std::cout << "m1list [4*N+5+"<< i <<"] in PV:" << std::endl;
        // m1list[4*N+5+i]->output();
    }
    temp8_1->sigma(temp9_1);
    temp8_1->mul(F[4] , temp9_1);
    F[4]->to_poly();
    // std::cout << "F in PV: " << std::endl;
    // F[4]->output();
    F[4]->polyarray[0] -= tempbetalist[1];

    // Polyq *temp8_2 = new Polyq(0);
    // Polyq *temp9_2 = new Polyq(0);
    // F[5] = new Polyq(0);
    // temp8_2->copy(m1list[0]);
    // temp8_2->sigma(temp9_2);
    // temp8_2->mul(F[5], temp9_2);
    // F[5]->to_poly();
    
    // std::cout << "F in PV: " << std::endl;
    // F[4]->output();
    
    //compute hvec
    polyvecq *hvec = new polyvecq(TAU);
    computeh(hvec, g, F, Gamma);
    // std::cout << "hvec in PV: " << std::endl;
    // hvec->polyarray[0]->output();
    //FourR step1

    // H3input
    
    unsigned char museed[TAU*POLYLEN];
    hvec->to_char(museed);
    keccak_state mu_state;
    shake256_init(&mu_state);
    shake256_absorb(&mu_state , Rseed, POLYLEN*(2+TAU+2*N1)+32+VECETALEN);
    shake256_absorb(&mu_state , Gammaseed , 256*8);
    shake256_absorb(&mu_state , museed , TAU*POLYLEN);

    polyvecq *mu = new polyvecq(&mu_state, 0, TAU , 0);

    
    //FourR step 2
    polymatq *Fymat = new polymatq(TAU+2 , M2);
    polyvecq *y = new polyvecq((M1_BATCH)*2 + (TAU+2)*2);
    
    constructFy_batch(Fymat ,Fg , F3);
    constructy_batch(y , Fymat , y1 , y2);
    
    //FourR step 3
    Polyq *D2list[N+4];
    for (int i = 0; i < N+4; i++)
    {
        D2list[i] = new Polyq(1);
    }
    
    constructD2_batch(D2list , mu , Gamma, c);
    // std::cout << "D2list in PV: " << std::endl;
    
    polyvecq *d1vec = new polyvecq((M1_BATCH)*2 + (TAU+2)*2);
    // d1vec->polyarray[(M1_BATCH)*2 + (TAU+2)*2-1]->output();

    constructd1_batch(d1vec , mu , Gamma , A , c);
    
    //FourR step 4
    
    //construct m
    polyvecq *m = new polyvecq((M1_BATCH)*2 + (TAU+2)*2);
    constructm_batch(m , m1 , g , y3 );
    

    // // test
    // Polyq *f = new Polyq(1);
    // Polyq *temp8_11 = new Polyq(1);
    // Polyq *temp8_12 = new Polyq(1);
    // Polyq *d0 = new Polyq(1);
    // mu->mul(f , hvec);
    // D2mul2(temp8_11 , D2list , m , m);
    // d1vec->mul(temp8_12 , m);
    // temp8_11->add(temp8_11 , temp8_12);
    // constructd2(d0 , mu , Gamma , h , hvec);
    // temp8_11->to_poly();
    // temp8_11->add(temp8_11 , d0);
    // temp8_11->output();

    // f->to_poly();
    // f->output();
    
    // compute g1
    Polyq tempg(1);
    Polyq *g1 = new Polyq(1);

    D2mul_batch(&tempg , D2list , m , y);

    g1->add(g1 , &tempg);
    D2mul_batch(&tempg , D2list,y , m);
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
    D2mul_batch(&tempg , D2list , y , y);
    v->add(v , &tempg);
    fg->mul(&tempg , y2);
    v->add(v , &tempg);
    v->to_poly();
    // std::cout << "v in PV: " << std::endl;

    //FifthR

    //H4input
    unsigned char cseed[POLYLEN*2];

    tpoly->to_char(cseed);
    v->to_char(cseed+POLYLEN);
    keccak_state c_state;
    shake256_init(&c_state);
    shake256_absorb(&c_state , Rseed, POLYLEN*(2+TAU+2*N1)+32+VECETALEN);
    shake256_absorb(&c_state , Gammaseed , 256*8);
    shake256_absorb(&c_state , museed , TAU*POLYLEN);
    shake256_absorb(&c_state , cseed , POLYLEN*2);
    Polyq *cpoly = new Polyq(0);
    Csamplepoly(&c_state , cpoly);

    // std::cout<<"seeds in PV:"<<std::endl;
    // printfBstr(Rseed,10);
    // printfBstr(Gammaseed,10);
    // printfBstr(museed,10);
    // printfBstr(cseed,10);

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
    // std::cout << "cpoly in PV: " << std::endl;
    // cpoly->output();

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
    // std::cout << z1->norm() << " " << z2->norm() << " " << vec_norm(z3vec , 256) << " " << std::endl;
    // std::cout << B1*B1 << " " << B2*B2 << " " << B3*B3 << " " << std::endl;
    // std::cout << (z1->norm() >B1*B1) << " " << (z2->norm() > B2*B2) << " " << (vec_norm(z3vec , 256) > B3*B3) << " " << std::endl;
    // if(z1->norm() > B11*B11 || z2->norm() > B22*B22 || vec_norm(z3vec , 256) > B33*B33){
    //     return 0;
    // }
    memcpy(pi , Rseed+VECETALEN+32 , POLYLEN*(2+TAU+2*N1));
    memcpy(pi+POLYLEN*(2+TAU+2*N1) , Gammaseed,256*8);
    memcpy(pi + POLYLEN*(2+TAU+2*N1) + 256*8 , museed, TAU*POLYLEN);

    
    memcpy(pi + POLYLEN*(2+TAU+2*N1) + 256*8 + TAU*POLYLEN, cseed , POLYLEN*2);

    // std::cout<<"z1 z2 in PV" <<std::endl;
    // z1->polyarray[0]->output();
    // z2->polyarray[0]->output();

    z1->to_char(pi + POLYLEN*(2+TAU+2*N1) + 256*8 + TAU*POLYLEN + POLYLEN*2); //z1 length: M1_BATCH*Polylen
    z2->to_char(pi + POLYLEN*(2+TAU+2*N1+M1_BATCH) + 256*8 + TAU*POLYLEN + POLYLEN*2); // z2 length: M2*polylen

    free(R);
    delete F3 ;
    delete Fg ;
    delete fg ;
    delete  E1 ;
    delete E2 ;
    delete cm1;
    delete temp8;
    delete temp8_1;
    delete temp9;
    delete temp9_1;
    free(m1vec);
    delete(w);
// POLYLEN*(4+TAU+2*N1) + 32+VECETALEN + 256*8 + TAU*POLYLEN + POLYLEN*2 + M1_BATCH + M2
    return 1;
}

int Verify_batch(unsigned char *crsseed, unsigned char *pi , polymatq *A, polymatq *c , polyvecq *t0 , polyvecq* t1 , polyvecq *h){
    keccak_state crsstate;
    shake256_init(&crsstate);
    shake256_absorb(&crsstate , crsseed,32);

    polymatq *E1 = new polymatq(&crsstate, 0, N1 ,M1_BATCH , 1);
    polymatq *E2 = new polymatq(&crsstate, 0, N1 ,M2 , 1);
    polymatq *Fg = new polymatq(&crsstate, 0, TAU ,M2, 1);
    polymatq *F3 = new polymatq(&crsstate, 0, 2 ,M2 , 1);
    polymatq *F4 = new polymatq(&crsstate, 0, 2 ,M2 , 1);
    polyvecq *fg = new polyvecq(&crsstate, 0, M2, 1);
    
    unsigned char Rseed[POLYLEN*(2+TAU+2*N1)+32+VECETALEN];
    memcpy(Rseed , crsseed , 32);

    //todo: insert statement x
    memset(Rseed+32 , 0 , VECETALEN);

    memcpy(Rseed+VECETALEN+32 , pi , POLYLEN*(2+TAU+2*N1));

    // unpack alpha1
    polyvecq *t3 = new polyvecq(2);
    //polyvecq *t4 = new polyvecq(2);
    polyvecq *tg = new polyvecq(TAU);
    polyvecq *tE = new polyvecq(N1);
    polyvecq *w = new polyvecq(N1);
    t3->from_char(Rseed+32+VECETALEN , 0);
    //t4->from_char(Rseed+32+VECETALEN+POLYLEN*2 , 0);
    tg->from_char(Rseed+32+VECETALEN+POLYLEN*2 , 0);
    tE->from_char(Rseed+32+VECETALEN+POLYLEN*(2+TAU) , 0);
    w->from_char(Rseed+32+VECETALEN+POLYLEN*(2+TAU+N1) , 0);

    // std::cout << "E1 in VF: " << std::endl;
    // E1->vecarray[0]->polyarray[0]->output();
    // std::cout << "E2 in VF: " << std::endl;
    // E2->vecarray[0]->polyarray[0]->output();


    keccak_state Rstate;
    shake256_init(&Rstate);
    shake256_absorb(&Rstate , Rseed, POLYLEN*(2+TAU+2*N1)+32+VECETALEN);
    
    int64_t *R = (int64_t*)malloc(sizeof(int64_t)*256*(M1_BATCH)*128);

    cbd(&Rstate , R , 1 , 256*128*(3*N+5));
    clock_t time1 = clock();
    
    //unpack alpha2

    int64_t z3vec[256];
    unsigned char Gammaseed[256*8];
    memcpy(Gammaseed ,pi+POLYLEN*(2+TAU+2*N1) ,256*8);
    char2vec(z3vec , Gammaseed , 256);

    keccak_state gamma_state;
    shake256_init(&gamma_state);
    shake256_absorb(&gamma_state , Rseed,POLYLEN*(2+TAU+2*N1)+32+VECETALEN);
    shake256_absorb(&gamma_state , Gammaseed , 256*8);
    
    polymatq *Gamma = new polymatq(TAU , N_BATCH+2*N+6);
    sample_gamma(&gamma_state , Gamma);


    // unpack alpha3
    unsigned char museed[TAU*POLYLEN];
    memcpy(museed, pi + POLYLEN*(2+TAU+2*N1) + 256*8 , TAU*POLYLEN);
    polyvecq *hvec = new polyvecq(TAU);
    hvec->from_char(museed , 0);
    keccak_state mu_state;
    shake256_init(&mu_state);
    shake256_absorb(&mu_state , Rseed, POLYLEN*(2+TAU+2*N1)+32+VECETALEN);
    shake256_absorb(&mu_state , Gammaseed , 256*8);
    shake256_absorb(&mu_state , museed , TAU*POLYLEN);

    polyvecq *mu = new polyvecq(&mu_state ,0, TAU , 0);

    //unpack alpha4
    unsigned char cseed[POLYLEN*2];
    memcpy(cseed , pi + POLYLEN*(2+TAU+2*N1) + 256*8 + TAU*POLYLEN , POLYLEN*2);
    Polyq *tpoly = new Polyq(0);
    Polyq *v = new Polyq(0);
    tpoly->from_char(cseed , 0);
    v->from_char(cseed+POLYLEN , 0);
    keccak_state c_state;
    shake256_init(&c_state);
    shake256_absorb(&c_state , Rseed, POLYLEN*(2+TAU+2*N1)+32+VECETALEN);
    shake256_absorb(&c_state , Gammaseed , 256*8);
    shake256_absorb(&c_state , museed , TAU*POLYLEN);
    shake256_absorb(&c_state , cseed , POLYLEN*2);
    Polyq *cpoly = new Polyq(0);
    Csamplepoly(&c_state , cpoly);


    // std::cout<<"seeds in VF:"<<std::endl;
    // printfBstr(Rseed,10);
    // printfBstr(Gammaseed,10);
    // printfBstr(museed,10);
    // printfBstr(cseed,10);

    //constructz
    polyvecq *z1 = new polyvecq(M1_BATCH);
    polyvecq *z2 = new polyvecq(M2);
    polyvecq *z = new polyvecq(2*(M1_BATCH)+2*(2+TAU));
    z1->from_char(pi + POLYLEN*(2+TAU+2*N1) + 256*8 + TAU*POLYLEN + POLYLEN*2 , 0);
    z2->from_char(pi + POLYLEN*(2+TAU+2*N1+M1_BATCH) + 256*8 + TAU*POLYLEN + POLYLEN*2 , 0);
    
    polymatq *Fymat = new polymatq(TAU+2 , M2);
    constructFy_batch(Fymat ,Fg , F3);
    constructz_batch(z , Fymat , tg , t3, cpoly , z1 , z2);

    // constructD
    Polyq *D2list[4+N];
    for (int i = 0; i < N+4; i++)
    {
        D2list[i] = new Polyq(1);
    }
    constructD2_batch(D2list , mu , Gamma, c);
    polyvecq *d1vec = new polyvecq(2*(M1_BATCH)+2*(2+TAU));
    constructd1_batch(d1vec , mu , Gamma , A, c);

    Polyq *d0 = new Polyq(1);
    constructd0_batch(d0 , mu , Gamma , h , hvec, t0, t1);
    //check norm
    // std::cout << z1->norm() << " " << z1->norm() << " " << vec_norm(z3vec , 256) << " " << std::endl;
    if(z1->norm() > B11*B11 || z2->norm() > B22*B22 || vec_norm(z3vec , 256) > B33*B33){
        std::cout << "wrong1" << std::endl;
        // return 0;
    }
    // std::cout<<"z1 z2 in VF" <<std::endl;
    // z1->polyarray[0]->output();
    // z2->polyarray[0]->output();
    // check E1z1 + E2z2 == w+ ctE

    polyvecq left(N1 , 1) , temp(N1 , 1) , right(N1 , 1);

    E1->right_mul(&temp , z1);
    left.add(&left, &temp); //E1*z1
    
    E2->right_mul(&temp , z2);
    left.add(&left , &temp);

    // std::cout << "condintion 2 left: E1z1+E2z2 in VF" << std::endl;
    // left.polyarray[0]->output(); //E2*z2

    w->to_ntt();

    //test w
    // std::cout << "w in VF: " << std::endl;
    // w->polyarray[0]->output();

    // std::cout << "tE in VF: " << std::endl;
    // tE->polyarray[0]->output();
    //  std::cout << "cpoly in VF: " << std::endl;
    //  cpoly->output();

    tE->mul_poly(&right , cpoly);
    right.add(&right , w);
    
    //test w+ ctE
    std::cout << "condintion 2 left: w+ctE in VF: " << std::endl;
    // right.to_poly();
    // right.polyarray[0]->output();

    if(left.equal(&right) == 0){
        std::cout << "wrong2" << std::endl;// return E1z1 + E2z2 == w+ ctE
        return 0;
    }

        //check h
    hvec->to_poly();
    // hvec->polyarray[0]->output();
    for(int i = 0 ; i < TAU ;i++){
        if(hvec->polyarray[i]->polyarray[0] !=0 ){
            std::cout << "wrong4" << std::endl;
            return 0;
        }
    }
    
    // check z^TD2z + cd1^Tz + c^2d0 + fg^Tz2= v + ct
    
    Polyq leftpoly(1) , temppoly(1);
    
    D2mul_batch(&leftpoly , D2list ,z ,  z);
    
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
    // if(leftpoly.equal(v) == 0){
    //     std::cout << "wrong3" << std::endl;
    //     return 0;
    // }
        
    delete F3 ;
    delete Fg ;
    delete fg ;
    delete E1 ;
    delete E2 ;
    delete w;
    free(R);
    return 1;
}
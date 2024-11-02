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
#include "protocol_batch.h"
#include "cpaenc.h"
#include "tools.h"
#include <ctime>
double max(double list[] , int len){
    double res = list[0];
    for(int i = 0 ; i < len ; i++){
        if(list[i] > res){
            res = list[i];
        }
    }
    return res;
}
double min(double list[] , int len){
    double res = list[0];
    for(int i = 0 ; i < len ; i++){
        if(list[i] < res){
            res = list[i];
        }
    }
    return res;
}
double average(double list[] , int len){
    double res = 0;
    for(int i = 0 ; i < len ; i++){
        res += list[i];
    }
    res /= len;
    return res;
}

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
    packsk(sk , s0, s1, k);

    A->trans();
       
    A->right_mul(t0 , s0);
    
    e0->to_ntt();
    t0->add(t0 , e0);

    A->right_mul(t1 , s1);
    
    e1->to_ntt();
    t1->add(t1 , e1);
    
    packpk(pk , seed1, t0 , t1);

}

void AT_CQ(unsigned char* pi1 ,unsigned char*crsseed, unsigned char* query, unsigned char* H1input ,\
unsigned char* xchar , unsigned char* pk){

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

    //unsigned char H1input[256+MD];
    memcpy(H1input , tag , MD/8);
    shake256_init(&state);
    shake256_absorb(&state , H1input , MD/8);
    
    polyvecq *ctemp = new polyvecq(&state , 0 ,N, 0);
//remove H(tag)
    // polyvecq *ctemp = new polyvecq(N,1);
    // std::cout<< "hash in IT:" <<std::endl;
    // ctemp->polyarray[0]->output();
// remove e2
    ctemp->add(ctemp , e2);
    ctemp->to_ntt();
    polyvecq *c = new polyvecq(N,1);
// remove Ax
    A->right_mul(c , x);
    c->add(c , ctemp);
    // std::cout << "ci in CQ: " << std::endl;
    // c->polyarray[0]->output();
    delete ctemp;
    c->to_char(query);//ntt
    
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
    // std::cout<< "x in cq:" << x->polyarray[0]->polyarray[1] <<std::endl;

    clock_t t2 = clock();
    //std::cout << "CQ " << (double)(t2-t1)/CLOCKS_PER_SEC << std::endl;

    while(1){
        if(Prove1(pi1 , crsseed , rseed , x , e2 , H1input , uv , pke , A , c) == 1){
            break;
        }
    }
    // std::cout << "pi1 in CQ:" << std::endl;
    // printfBstr(pi1,10);
    clock_t t3 = clock();
    //std::cout <<"prove1 " << (double)(t3-t2)/CLOCKS_PER_SEC << std::endl;
    std::cout << "CQ:done" << std::endl;
}

void AT_IT(unsigned char* pi2 , unsigned char* crsseed2 , unsigned char* resp ,unsigned char* pk,\
unsigned char *sk, unsigned char** query , unsigned char** crsseed , unsigned char** pi1, double* pvtime){
    clock_t time1 = clock();
    polyvecq *s0 = new polyvecq(N);
    polyvecq *s1 = new polyvecq(N);
    polyvecq *s2 = new polyvecq(N);
    polyvecq *e0 = new polyvecq(N);
    polyvecq *e1 = new polyvecq(N);
    polymatq *c = new polymatq(N_BATCH,N);
    polyvecq *t0 = new polyvecq(N);
    polyvecq *t1 = new polyvecq(N);
    unsigned char seed1[32];
    unsigned char seedA[32];
    unsigned char H2input[VECLEN*N_BATCH+32];
    int bit[N_BATCH];
    memcpy(seedA , pk , 32);
    memcpy(seed1 , pk+32 , 32);
    keccak_state state;

    keccak_state Astate;
    shake256_init(&Astate);
    shake256_absorb(&Astate , seedA , 32);
    
    polymatq *A = new polymatq(&Astate ,0, N , N , 1);

    t0->from_char(pk+32 , 1);
    t1->from_char(pk+32+VECLEN , 1);
    unpacksk(sk, s0, s1, H2input);
    //Polyq* temp;
    polyvecq* e3 = new polyvecq(N_BATCH);
    polyvecq* h = new polyvecq(N_BATCH);

    memcpy(H2input+32 , query , VECLEN);

    shake256_init(&state);
    shake256_absorb(&state , H2input , VECLEN*N_BATCH+32);
    e3 = new polyvecq(&state , ETAE , N_BATCH, 0);

    for(int i = 0 ; i < N_BATCH; i++){
        c->vecarray[i]->from_char(query[i] , 1); //ntt
        bit[i] = rand() % 2;
        if(bit[i] == 1){
            s1->mul(h->polyarray[i] , c->vecarray[i]);
            // std::cout<< "1" <<std::endl;
        }
        else{
            s0->mul(h->polyarray[i] , c->vecarray[i]);
            // std::cout<< "0" <<std::endl;
        }
        h->polyarray[i]->to_poly();
//remove error term
        h->polyarray[i]->add(h->polyarray[i],e3->polyarray[i]);

        // std::cout<< "hi in IT:" <<std::endl;
        // h->polyarray[i]->output();

        h->polyarray[i]->to_char(resp+i*POLYLEN);
        // std::cout<< "bit in IT:" << bit[i] <<std::endl;

        // std::cout << " pi1 in IT: "<<std::endl;
        // c->vecarray[i]->polyarray[0]->output();
        // printfBstr(pi1[i],10);
        std::cout << "verify1 "<< Verify1(crsseed[i], pi1[i] , A , c->vecarray[i]) << std::endl;
    }
    clock_t t2 = clock();
    
    //output time
    // std::cout <<"IT " << (double)(t2-time1)/CLOCKS_PER_SEC << std::endl;
    
    setup(crsseed2);
    
    int loop = 0;
    while(!loop){
         loop = Prove_batch(pi2 , crsseed2 , bit , s0, s1 , e0 , e1 , e3 , h , A , t0, t1 , c);
    }
    clock_t t4 = clock();
    if(pvtime != NULL){
        pvtime[0] = (double)(t4-t2)/CLOCKS_PER_SEC;
    }

    // std::cout <<"prove_batch " << (double)(t4-t2)/CLOCKS_PER_SEC << std::endl;
}

void AT_CF(unsigned char** msg , unsigned char* pk , unsigned char*resp ,unsigned char**H1input, unsigned char** xchar ,\
    unsigned char** query , unsigned char* pi2 , unsigned char* crsseed2, double* vftime){
    
    clock_t time1 = clock();
    polyvecq *t0 = new polyvecq(N);
    polyvecq *t1 = new polyvecq(N);  
    polymatq *x = new polymatq(N_BATCH, N);
    polyvecq* h = new polyvecq(N_BATCH,0);
    polyvecq* v0 = new polyvecq(N_BATCH,0);
    polyvecq* v1 = new polyvecq(N_BATCH,0);
    polyvecq* temp1 = new polyvecq(N_BATCH,1), *temp2 = new polyvecq(N_BATCH,1); //ntt
    unsigned char v0char[N_BATCH][COMPOLYLEN] , v1char[N_BATCH][COMPOLYLEN];
    
    t0->from_char(pk+32 , 1);
    t1->from_char(pk+32+VECLEN , 1);

//extract h1, ..., hN
    for (int i = 0; i < N_BATCH; i++){  
        h->polyarray[i]->from_char(resp+i*POLYLEN, 0);
        // std::cout<< "hi in cf:" <<std::endl;
        // h->polyarray[i]->output();
    }

//recover x
    for (int i = 0; i < N_BATCH; i++){
        // std::cout<< "x in cf:" <<std::endl;
        x->vecarray[i]-> from_char(xchar[i],0,ETAR);
        t0->mul(temp1->polyarray[i] , x->vecarray[i]);
        
        temp1->polyarray[i]->to_poly();
        // std::cout<< "<t0,x> in cf:" <<std::endl;
        // temp1->polyarray[i]->output();

        h->polyarray[i]->sub(v0->polyarray[i] , temp1->polyarray[i]);
        // std::cout<< "v0 in cf:" <<std::endl;
        // v0->polyarray[i]->output();
        Compress(v0char[i] , v0->polyarray[i]);
        //  std::cout<< "v0char in cf:" <<std::endl;
        //  printfBstr(v0char[i],32);

        t1->mul(temp2->polyarray[i] , x->vecarray[i]);
        temp2->polyarray[i]->to_poly();
        // std::cout<< "<t1,x> in cf:" <<std::endl;
        // temp1->polyarray[i]->output();
        h->polyarray[i]->sub(v1->polyarray[i] , temp2->polyarray[i]);
        // std::cout<< "v1 in cf:" <<std::endl;
        // v1->polyarray[i]->output();
        Compress(v1char[i] , v1->polyarray[i]);
        //  std::cout<< "v1char in cf:" <<std::endl;
        //  printfBstr(v1char[i],32);
    }

//store N_BATCH tuples (v0,v1,tag)  
    for (int i = 0; i < N_BATCH; i++){
        memcpy(msg[i] , H1input[i] ,32);
        memcpy(msg[i]+32 , v0char[i] , COMPOLYLEN);
        memcpy(msg[i]+32+COMPOLYLEN , v1char[i] , COMPOLYLEN);
        
        // std::cout <<"msg in CF " << std::endl;
        // printfBstr(msg[i],10);
    }
//recover c
    polymatq *c = new polymatq(N_BATCH, N);
    for (int i = 0; i < N_BATCH; i++){
        c->vecarray[i]->from_char(query[i] , 1);
    }

    keccak_state Astate;
    shake256_init(&Astate);
    shake256_absorb(&Astate , pk , 32);

    polymatq *A = new polymatq(&Astate ,0, N , N , 1);

    clock_t t2 = clock();
    // std::cout <<"CF " << (double)(t2-time1)/CLOCKS_PER_SEC << std::endl;
    std::cout << "verify_batch " << Verify_batch(crsseed2 , pi2 , A , c , t0 , t1 , h) << std::endl;
    clock_t t3 = clock();
    // std::cout <<"verify_batch " << (double)(t3-t2)/CLOCKS_PER_SEC << std::endl;
    if(vftime != NULL){
        vftime[0] = (double)(t3-t2)/CLOCKS_PER_SEC;
    }
}
int AT_Rb(unsigned char* pk , unsigned char* sk , unsigned char* msg){

    clock_t t1 = clock();
    keccak_state state;
    int bit = 0;
    unsigned char H1input[32];
    unsigned char res[COMPOLYLEN];
    memcpy(H1input , msg , 32);
    // printfBstr(H1input,32);
    
    shake256_init(&state);
    shake256_absorb(&state , H1input , 32);
    polyvecq *h = new polyvecq(&state , 0 ,N, 0);
    // std::cout<< "hash in rb:" <<std::endl;
    // h->polyarray[0]->output();
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
    
    // std::cout<< "v0 in rb:" <<std::endl;
    // // printfBstr(msg+32,32);
    // v0_->output();
    // std::cout<< "v1 in rb:" <<std::endl;
    // // printfBstr(msg+32+COMPOLYLEN,32);
    // v1_->output();


    s0->mul(u0 , h);
    s1->mul(u1 , h);
    
    u0->to_poly();
    u1->to_poly();
    // std::cout<< "u0 in rb:" <<std::endl;
    // u0->output();
    // std::cout<< "u1 in rb:" <<std::endl;
    // u1->output();

    u0->sub(temp , v0_);
    // std::cout<< "u0-v0 in rb:" <<std::endl;
    // temp->output();
    Compress(res , temp);
    // printfBstr(res , COMPOLYLEN);
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
    // std::cout<< "u1-v1 in rb:" <<std::endl;
    // temp->output();
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

int main(){
    unsigned char pk[PKLEN] , sk[SKLEN];
    unsigned char *crsseed[N_BATCH], *H1input[N_BATCH] , *xchar[N_BATCH];
    unsigned char resp[VECLEN*N_BATCH] , crsseed2[32];
    unsigned char* query[N_BATCH], *pi1[N_BATCH];
    unsigned char* msg[N_BATCH];
    int bit[N_BATCH];
    unsigned char* pib = (unsigned char*)malloc(PIBLEN);
    keygen(pk , sk);
    int n = 5;
    double ITtime[n], pvtime[n], vftime[n];
    time_t t1 , t2;

    for(int j = 0 ; j < n ; j++){
        for(int i=0;i<N_BATCH;i++){
            msg[i] = (unsigned char*)malloc(VECLEN+32+MD/8);
            xchar[i] = (unsigned char*)malloc(VECETALEN);
            crsseed[i] = (unsigned char*)malloc(32);
            H1input[i] = (unsigned char*)malloc(32+MD/8);
        }
        for(int i = 0; i < N_BATCH ; i++){
            pi1[i] = (unsigned char*)malloc(PI1LEN);
            query[i] = (unsigned char*)malloc(VECLEN);
            AT_CQ(pi1[i] ,crsseed[i], query[i] , H1input[i] , xchar[i] , pk);
        }
        t1 = clock();
        AT_IT(pib , crsseed2 , resp , pk , sk , query , crsseed , pi1, pvtime+j);
        t2 = clock();
        ITtime[j] = (double)(t2-t1)/CLOCKS_PER_SEC;
        AT_CF(msg , pk , resp , H1input , xchar , query , pib , crsseed2, vftime+j);
        for(int i=0; i < N_BATCH; i++){
            bit[i] = AT_Rb(pk , sk , msg[i]);
            std::cout << " Readbit: " << bit[i] << std::endl;
        }
        for(int i=0;i<N_BATCH;i++){
            free(msg[i]);
            free(xchar[i]);
            free(crsseed[i]);
            free(H1input[i]);
            free(pi1[i]);
            free(query[i]);
        }
        
    }

    std::cout << max(ITtime , n) <<" "<< min(ITtime , n) << " " << average(ITtime , n) << std::endl;
    std::cout << max(pvtime , n) <<" "<< min(pvtime , n) << " " << average(pvtime , n) << std::endl;
    std::cout << max(vftime , n) <<" "<< min(vftime , n) << " " << average(vftime , n) << std::endl;
    // std::cout << " multi-protocol: done " << std::endl;
    
    return 0;
}
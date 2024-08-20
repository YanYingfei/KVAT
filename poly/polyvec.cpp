#include <stdlib.h>
#include "polymat.h"

template <typename T>
polyvec<T>::polyvec(T *vec[] , int k):k(k){
    polyarray = (T**)malloc(k*sizeof(T*));
    for(int i = 0 ; i < k ;i++){
        polyarray[i] = vec[i];
    }
}


template <typename T>
polyvec<T>::polyvec(int k , int nttflag):k(k){
    polyarray = (T**)malloc(k*sizeof(T*));
    for(int i = 0 ; i < k ;i++){
        polyarray[i] = new T(nttflag);
    }
}

template <typename T>
polyvec<T>::polyvec(keccak_state *state, int eta, int k, int nttflag):k(k){
    polyarray = (T**)malloc(k*sizeof(T*));
    for(int i = 0 ; i < k ; i++){
        polyarray[i] = new T(state , eta , nttflag);
    }
}

template <typename T>
polyvec<T>::polyvec(unsigned char *seed, int seedlen, int eta, int k, int nttflag):k(k){
    keccak_state state;
    shake256_init(&state);
    shake256_absorb(&state , seed , seedlen);
    polyarray = (T**)malloc(k*sizeof(T*));
    for(int i = 0 ; i < k ; i++){
        polyarray[i] = new T(&state , eta , nttflag);
    }
}
template <typename T>
polyvec<T>::polyvec(double sigma, int k, int nttflag):k(k){
    polyarray = (T**)malloc(k*sizeof(T*));
    for(int i = 0 ; i < k ; i++){
        polyarray[i] = new T(sigma, nttflag);
    }
}

template <typename T>
void polyvec<T>::reset(int nttflag){
    for(int i = 0 ; i < this->k ; i++){
        this->polyarray[i]->reset(nttflag);
    }
}

template <typename T>
void polyvec<T>::mul(T *res, polyvec *rvalue)
{
    assert(this->k == rvalue->k);
    T *mid = new T(1);
    res->reset(1);
    for(int i = 0 ; i < k;i++){
        this->polyarray[i]->mul(mid , rvalue->polyarray[i]);
        res->add(res , mid);
    }
    delete mid;
}
template <typename T>
void polyvec<T>::add(polyvec *res, polyvec *rvalue){
    assert(this->k == rvalue->k);
    for(int i = 0 ; i < k ; i++){
        this->polyarray[i]->add(res->polyarray[i] , rvalue->polyarray[i]);
    }
}

template <typename T>
void polyvec<T>::sub(polyvec *res, polyvec *rvalue){
    assert(this->k == rvalue->k);
    for(int i = 0 ; i < k ; i++){
        this->polyarray[i]->sub(res->polyarray[i] , rvalue->polyarray[i]);
    }
}
template <typename T>
int polyvec<T>::equal(polyvec *rvalue){
    assert(this->k == rvalue->k);
    for(int i = 0 ; i <this->k ; i++){
        if(this->polyarray[i]->equal(rvalue->polyarray[i]) != 1){
            return 0;
        }
    }
    return 1;
}
template <typename T>
void polyvec<T>::mul_poly(polyvec *res, T *rvalue)
{
    assert(this->k == res->k);
    for(int i = 0 ; i < this->k ; i++){
        this->polyarray[i]->mul(res->polyarray[i] , rvalue);
    }
}
template <typename T>
void polyvec<T>::to_poly()
{
    for(int i = 0 ; i < k ; i++){
        this->polyarray[i]->to_poly();
    }
}
template <typename T>
void polyvec<T>::to_ntt(){
    for(int i = 0 ; i < k ; i++){
        this->polyarray[i]->to_ntt();
    }
}

template <typename T>
void polyvec<T>::output(){
    for(int i = 0 ; i < k ; i++){
        this->polyarray[i]->output();
    }
}

template <typename T>
void polyvec<T>::to_char(unsigned char *res , int eta){
    int polylen;
    if(eta == 0){
        polylen = this->polyarray[0]->len;
    }
    else{
        assert(eta == 7);
        polylen = POLYETALEN;
    }
    for(int i = 0 ; i < k ; i++){
        this->polyarray[i]->to_char(res+polylen*i , eta);
    }
}
template <typename T>
void polyvec<T>::from_char(unsigned char *vecchar , int nttflag ,int eta){
    int polylen;
    if(eta == 0){
        polylen = this->polyarray[0]->len;
    }
    else{
        assert(eta == 7);
        polylen = POLYETALEN;
    }
    for(int i = 0 ; i < k ; i++){
        this->polyarray[i]->from_char(vecchar+polylen*i , nttflag , eta);
    }
}
template <typename T>
void polyvec<T>::sigma(polyvec *res){
    for(int i = 0 ; i < this->k ; i++){
        this->polyarray[i]->sigma(res->polyarray[i]);
    }
}
template <typename T>
int64_t polyvec<T>::norm()
{
    int64_t res = 0;
    for(int i = 0 ; i < k; i++){
        res += polyarray[i]->norm();
    }
    return res;
}

template <typename T>
polyvec<T>::~polyvec()
{
    for(int i = 0;i<k;i++){
        delete polyarray[i];
    }
    free(polyarray);
}
template class polyvec<Polyq>;
template class polyvec<Polyp>;
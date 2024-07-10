#include <stdlib.h>
#include <iostream>
#include "polymat.h"
template <typename T>
polymat<T>::polymat(polyvec<T> *vecarray[] , int k , int l):k(k),l(l){
    this->vecarray = (polyvec<T>**)malloc(k*sizeof(polyvec<T>**));
    for(int i = 0;i<k;i++){
        this->vecarray[i] = vecarray[i];
    }
}
template <typename T>
polymat<T>::polymat(int k , int l):k(k),l(l){
    this->vecarray = (polyvec<T>**)malloc(k*sizeof(polyvec<T>**));
    for(int i = 0;i<k;i++){
        this->vecarray[i] = new polyvec<T>(l);
    }
}
template <typename T>
void polymat<T>::right_mul(polyvec<T> *res, polyvec<T> *rvalue){
    assert(this->l == rvalue->k);
    assert(this->k == res->k);
    for(int i = 0 ; i < k ;i++){
        this->vecarray[i]->mul(res->polyarray[i] , rvalue);
    }
}

template <typename T>
void polymat<T>::trans(){
    T *temp;
    //polyvec** newmat;
    if(k == l){
        for(int i = 0; i < k-1;i++){
            for(int j = i+1; j < k ; j++){
                temp = this->vecarray[i]->polyarray[j];
                this->vecarray[i]->polyarray[j] = this->vecarray[j]->polyarray[i];
                this->vecarray[j]->polyarray[i] = temp;
            }
        }
    }
    else{
        assert(1==0);
        //newmat = (polyvec**)malloc(k*sizeof(polyvec**));

    }
}
template <typename T>
polymat<T>::polymat(keccak_state *state, int eta, int k,int l, int nttflag):k(k),l(l){
    this->vecarray = (polyvec<T>**)malloc(k*sizeof(polyvec<T>**));
    for(int i = 0;i<k;i++){
        this->vecarray[i] = new polyvec<T>(state , eta , l , nttflag);
    }
}
template <typename T>
polymat<T>::polymat(unsigned char *seed, int seedlen, int eta, int k, int l, int nttflag):k(k),l(l){
    keccak_state state;
    shake256_init(&state);
    shake256_absorb(&state , seed , seedlen);
    this->vecarray = (polyvec<T>**)malloc(k*sizeof(polyvec<T>**));
    for(int i = 0;i<k;i++){
        this->vecarray[i] = new polyvec<T>(&state , eta , l , nttflag);
    }
}
template <typename T>
void polymat<T>::to_ntt()
{
    for(int i = 0;i<k;i++){
        vecarray[i]->to_ntt();
    }
}

template <typename T>
void polymat<T>::to_poly(){
    for(int i = 0;i<k;i++){
        vecarray[i]->to_poly();
    }
}

template <typename T>
inline polymat<T>::~polymat(){
    for(int i = 0;i<k;i++){
        delete vecarray[i];
    }
    free(vecarray);
}
template class polymat<Polyp>;
template class polymat<Polyq>;
#ifndef POLYVEC_H
#define POLYVEC_H
#include <assert.h>
#include "poly.h"
#include "fips202.h"

template <typename T>
class polyvec{
public:
    T** polyarray;
    int k;
    polyvec(T* vec[] ,int k);
    polyvec(int k , int nttflag = 0);
    polyvec(keccak_state* state ,int eta, int k , int nttflag);
    polyvec(unsigned char* seed , int seedlen , int eta , int k , int nttflag);
    polyvec(double sigma , int k , int nttflag);
    void reset(int nttflag);
    void mul(T*res , polyvec* rvalue);
    void add(polyvec *res , polyvec *rvalue);
    void sub(polyvec *res , polyvec *rvalue);
    int equal(polyvec*rvalue);
    void mul_poly(polyvec*res , T *rvalue);
    void output();
    void to_poly();
    void to_ntt();
    void to_char(unsigned char* res, int eta = 0);
    void from_char(unsigned char* vecchar, int nttflag , int eta = 0);
    void sigma(polyvec *res);
    int64_t norm();
    ~polyvec();
};
template <typename T>
class polymat{
public:
    polyvec<T> **vecarray;
    int k , l;  // k * l , row:k ,col:l
    polymat(polyvec<T> *vecarray[] , int k , int l);
    polymat(int k , int l);
    polymat(keccak_state* state ,int eta, int k ,int l, int nttflag);
    polymat(unsigned char* seed , int seedlen , int eta , int k ,int l, int nttflag);
    void right_mul(polyvec<T>*res , polyvec<T> *rvalue);
    void trans();
    void to_ntt();
    void to_poly();
    ~polymat();
};

#endif


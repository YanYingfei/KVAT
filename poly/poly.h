#ifndef POLY_H
#define POLY_H
#include "params.h"
#include <stdint.h>
#include <cstring>
#include "fips202.h"
#include "sample.h"



template <typename T>
class Poly{
public:
    T polyarray[128];
    T* ntttable;
    T q, q2inv;
    int nttflag;
    Poly(T *array , int nttflag);
    Poly(int nttflag);
    void set(T *array , int nttflag);
    void reset(int nttflag);
    void add(Poly *res , Poly *b);
    void sub(Poly *res , Poly *b);
    void mul(Poly *res , Poly *b);
    int equal(Poly *b);
    void mul_num(Poly *res , T b);
    void to_ntt();
    void to_poly();
    void sigma(Poly *res);
    void copy(Poly *target);
    
    T mod(T a,int flag = 1);
    T mulmod(T a , T b);
    void output(int flag = 0);
    int64_t norm();
    
};

class Polyq:public Poly<int64_t>{
public:
    int len = POLYLEN;
    Polyq(int64_t *array , int nttflag);
    Polyq(int nttflag);
    Polyq(unsigned char*seed, int seedlen, int eta,int nttflag);
    Polyq(keccak_state *state ,int eta, int nttflag);
    Polyq(double sigma , int nttflag);
    void from_char(unsigned char *polychar , int nttflag , int eta = 0);
    void to_char(unsigned char*res ,int eta = 0);
};

class Polyp:public Poly<int16_t>{
public:
    int len = POLYPLEN;
    Polyp(int16_t *array , int nttflag);
    Polyp(int nttflag);
    Polyp(unsigned char*seed, int seedlen, int eta,int nttflag);
    Polyp(keccak_state *state ,int eta, int nttflag);
    Polyp(double sigma , int nttflag);
    void from_char(unsigned char *polychar , int nttflag , int eta = 0);
    void to_char(unsigned char*res ,int eta = 0);
};

#endif
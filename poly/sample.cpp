#include "sample.h"
#include <random>
void uniformq(keccak_state *state , int64_t *res , int len){
    unsigned char buf[136];
    int index = 0 , bufindex = 0;
    int64_t a, b;
    shake256_squeezeblocks(buf , 1 , state);
    while(index < len){
        a = 0;
        b = 0;
        for(int j = 0 ; j < 4; j++){
            a <<= 8;
            a |= buf[bufindex];
            bufindex += 1;
        }
        a <<= 4;
        a |= buf[bufindex]>>4;
        b = buf[bufindex]&0xf;
        bufindex += 1;
        for(int j = 5 ; j < 9; j++){
            b <<= 8;
            b |= buf[bufindex];
            bufindex += 1;
        }

        if(a < Q && index < len){
            res[index] = a;
            index += 1;
        }
        if(b < Q && index < len){
            res[index] = b;
            index += 1;
        }
        if(bufindex == 135){
            shake256_squeezeblocks(buf , 1 , state);
            bufindex = 0;
        }
    }
}
void uniformp(keccak_state *state , int16_t *res){
    unsigned char buf[136];
    int index = 0 , bufindex = 0;
    int16_t a, b , q = 3329;
    shake256_squeezeblocks(buf , 1 , state);
    while(index < 128){
        a = ((buf[bufindex])<<4) | (buf[bufindex+1]>>4);
        b = ((buf[bufindex+1]&0xf)<<8) | buf[bufindex+2];
        if(a < q && index < 128){
            res[index] = a;
            index += 1;
        }
        if(b < q && index < 128){
            res[index] = b;
            index += 1;
        }
        bufindex += 3;
        if(bufindex == 135){
            shake256_squeezeblocks(buf , 1 , state);
            bufindex = 0;
        }
    }
}

template <typename T>
void cbd(keccak_state*state , T *res , int eta , int len){
    unsigned char buf[136];
    int index = 0 , bufindex = 0 , bitindex = 0;
    T a;
    shake256_squeezeblocks(buf , 1 , state);
    for(int i = 0 ; i < len ; i++){
        a = 0;
        for(int j = 0; j < eta ;j++){
            a += (buf[bufindex]>>bitindex)&1;
            a -= (buf[bufindex]>>(bitindex+1))&1;
            bitindex += 2;
            if(bitindex == 8){
                bufindex++;
                if(bufindex == 136){
                    shake256_squeezeblocks(buf , 1 , state);
                    bufindex = 0;
                }
                bitindex = 0;
            }
        }
        res[i] = a;
    }
}
template <typename T>
void Gaussian(T *res , double sigma){
    std::random_device rd;
    std::mt19937 generator(rd());
    std::normal_distribution<double> dis(0,sigma);
    for(int i = 0 ; i < 128 ; i++){
        res[i] = (int)dis(generator);
    }
}

template void cbd<int16_t>(keccak_state* , int16_t*, int,int);
template void cbd<int64_t>(keccak_state* , int64_t*, int,int);
template void Gaussian<int64_t>(int64_t *res , double sigma);
template void Gaussian<int16_t>(int16_t *res , double sigma);

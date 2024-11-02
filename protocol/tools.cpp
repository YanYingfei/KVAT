#include "tools.h"
#include <random>
#include "../poly/sample.h"
#include <iostream>
void randseed(unsigned char* seed , int seedlen){
    srand(clock());
    for(int i = 0 ; i < seedlen ; i++){
        seed[i] = rand();
    }
}
void matrix_mul(int64_t *res , int64_t *matrix , int64_t *vector , int k , int l){
    int64_t temp;
    for(int i = 0; i < k ;i++){
        temp = 0;
        for(int j = 0; j < l ; j++){
            temp += matrix[i*l+j]*vector[j];
        }
        res[i] = temp;
    }
}
void vec_add(int64_t *a , int64_t *b , int len){
    for(int i = 0 ; i < len ; i++){
        a[i] += b[i];
    }
}

void vec_mul_si(int64_t *a , int64_t x , int len){
    for(int i = 0 ; i < len ; i++){
        a[i] *=x;
    }
}
int64_t vec_mul(int64_t *a , int64_t *b , int len){
    int64_t res = 0;
    for(int i = 0 ; i < len ; i++){
        res +=a[i]*b[i];
    }
    return res;
}

int64_t vec_norm(int64_t *vec , int veclen){
    return vec_mul(vec , vec, veclen);
}

int Rej0(int64_t *z , int64_t *v , double s , double GAMMA){
    double a;
    int64_t b = 2*(rand()%2)-1;
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> dis(0.0, 1.0);
    a = dis(gen);
    double M_ = exp(1/2/GAMMA/GAMMA);
    M_ *= exp(-vec_norm(v ,256)/2/s/s) * cosh(b*vec_mul(z , v , 256)/s/s);

    M_ = 1/M_;
    if(a > M_){
        return 1;
    }
    return 0;

}

int Rej1(polyvecq *z, polyvecq *v, double s , double GAMMA){
    int64_t *zvec = (int64_t*)malloc(sizeof(int64_t)*z->k*128);
    int64_t *vvec = (int64_t*)malloc(sizeof(int64_t)*z->k*128);
    polyvec2vec(zvec , z);
    polyvec2vec(vvec , v);
    
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> dis(0.0, 1.0);
    double a = dis(gen);
    double M_ = exp(14/GAMMA+1/2/GAMMA/GAMMA);
    int res;
    double temp = (vec_norm(vvec , z->k*128) - 2* vec_mul(zvec , vvec , z->k*128)) / (2*s*s);

    M_ = 1/M_;
    M_ *= exp(temp);
    if(a > M_){
        res = 1;
    }
    else{
        res = 0;
    }
    free(zvec);
    free(vvec);
    return res;
}

void sample_gamma(keccak_state *state, polymatq *res){
    int64_t temp[res->k*res->l];
    uniformq(state , temp , res->k*res->l);
    for(int i = 0 ; i < res->k ; i++){
        for(int j = 0 ; j < res->l ; j++){
            res->vecarray[i]->polyarray[j]->polyarray[0] = temp[i*res->l+j];
        }
    }
    
}

void vec2char(unsigned char *res, int64_t *vec , int len){
    for(int i = 0 ; i < len ;i++){
        for(int j = 0 ; j < 8 ; j++){
            res[i*8+j] = ((vec[i])>>((7-j)*8))&0xff;
            // std::cout << vec[i] <<" " <<  (((vec[i])>>((7-j)*8))&0xff) << std::endl;
        }
        
        
    }
}
void char2vec(int64_t *res, unsigned char *vecchar, int len){
    for(int i = 0 ; i < len ; i++){
        res[i] = 0;
        for(int j = 0 ; j < 8 ; j++){
            res[i] += (int64_t)vecchar[i*8+j] << ((7-j)*8);
        }
        // if(res[i] < 0){
        //     res[i] += 1;
        // }
    }
}
void polyvec2vec(int64_t *res, polyvecq *y)
{
    y->to_poly();
    for(int i = 0 ; i < y->k ;i++){
        for(int j = 0 ; j < 128 ; j++){
        res[i*128+j] = y->polyarray[0]->mod(y->polyarray[i]->polyarray[j] , 0);
        }
    }
}
void geta(Polyq *ai , int norm){
    int a,b,c,d;
    int a2,b2,c2;
    for(a = 0 ; a*a<=norm;a++){
        a2 = a*a;
        for(b=a;a2+b*b<=norm;b++){
            b2 = a2+b*b;
            for(c=b;b2+c*c<=norm;c++){
                c2 = b2+c*c;
                d = (int)sqrt(norm - c2);
                if(c2+d*d == norm){
                    ai->polyarray[0] = a;
                    ai->polyarray[1] = b;
                    ai->polyarray[2] = c;
                    ai->polyarray[3] = d;
                }
            }
        }
    }
}

int Cjudge(Polyq* a){
    Polyq temp(0);
    a->sigma(&temp);
    int flag = 1;
    for(int i = 0 ; i < 128 ; i++){
        if(temp.polyarray[i] != a->polyarray[i]){
            flag = 0;
            break;
        }
    }
    return flag;
}

void Csamplepoly(keccak_state*state , Polyq *res){
    while(1){
        cbd(state , res->polyarray , 2 , 64);
        for(int i = 1; i < 64;i++){
            res->polyarray[128-i] = -res->polyarray[i];
        }
        if(Cjudge(res)){
            break;
        }
    }
}

void printfBstr(unsigned char *s , int l){
    for(int i = 0; i < l ; i++){
        printf("%02X" , s[i]);
    }
    printf("\n");
}
#include "poly.h"

#include <iostream>
#include <cassert>


template <typename T>
Poly<T>::Poly(T *array, int nttflag):nttflag(nttflag){
    std::memcpy(this->polyarray , array , 128*sizeof(T));
}
template <typename T>
Poly<T>::Poly(int nttflag):nttflag(nttflag){
    std::memset(this->polyarray , 0 , 128*sizeof(T));
}

template <typename T>
void Poly<T>::set(T*array , int nttflag){
    std::memcpy(this->polyarray , array , 128*sizeof(T));
    this->nttflag = nttflag;
}
template <typename T>
void Poly<T>::reset(int nttflag){
    std::memset(this->polyarray , 0 , 128*sizeof(T));
    this->nttflag = nttflag;
}
template <typename T>
void Poly<T>::add(Poly *res, Poly *b){
    assert(this->nttflag == b->nttflag);
    for(int i = 0 ; i < 128 ; i++){
        res->nttflag = this->nttflag;
        res->polyarray[i] = mod(this->polyarray[i] + b->polyarray[i]);
    }
}
template <typename T>
void Poly<T>::sub(Poly *res , Poly *b){
    assert(this->nttflag == b->nttflag);
    for(int i = 0 ; i < 128 ; i++){
        res->nttflag = this->nttflag;
        res->polyarray[i] = mod(this->polyarray[i] - b->polyarray[i]);
    }
}
template <typename T>
void Poly<T>::mul(Poly *res, Poly *b){
    if(this->nttflag == 0){
        this->to_ntt();
    }
    if(b->nttflag == 0){
        b->to_ntt();
    }
    res->nttflag = 1;
    for(int i = 0 ; i < 128; i++){
        res->polyarray[i] = mulmod(this->polyarray[i] , b->polyarray[i]);
    }
}
template <typename T>
int Poly<T>::equal(Poly *b){
    for(int i = 0 ; i < 128;i++){
        if(this->polyarray[i] != b->polyarray[i]){
            return 0;
        }
    }
    return 1;
}
template <typename T>
void Poly<T>::mul_num(Poly *res, T b)
{
    this->to_poly();
    res->nttflag = 0;
    for(int i = 0 ; i < 128 ; i++){
        res->polyarray[i] = mulmod(this->polyarray[i] , b);
    }
}
template <typename T>
void Poly<T>::to_ntt()
{
    unsigned int len, start, j, k;
    T zeta;
    T t;
    if(nttflag){
        return;
    }
    nttflag = 1;
    
    k = 0;

    for(len = 64; len > 0; len >>= 1) {
        for(start = 0; start < 128; start += 2*len){
            zeta = ntttable[++k];
            for(j = start; j < start + len; ++j) {
                t = mulmod(zeta , polyarray[j + len]);
                polyarray[j + len] = mod(polyarray[j] - t);
                polyarray[j] = mod(polyarray[j] + t);
            }
        }
    }
}
template <typename T>
void Poly<T>::to_poly(){
    unsigned int start, len, j, k;
    T t;
    T zeta;
    if(nttflag == 0){
        return ;
    }
    nttflag = 0;
    k = 128;

    for(len = 1; len < 128; len <<= 1) {
        for(start = 0; start < 128; start += 2* len) {
            zeta = -ntttable[--k];
            for(j = start; j < start + len; ++j) {
                t = polyarray[j];
                polyarray[j] = mod(t + polyarray[j + len]);
                polyarray[j + len] = mod(t - polyarray[j + len]);
                polyarray[j + len] = mulmod(zeta , polyarray[j + len]);
            }
        }
    }
    for(int i = 0 ; i < 128; i++){
        polyarray[i] = mulmod(q2inv,polyarray[i]);
    }
}
template <typename T>
void Poly<T>::sigma(Poly *res){
    res->set(this->polyarray , this->nttflag);
    if(res->nttflag == 1){
        res->to_poly();
    }
    T a;
    for(int i = 1; i <= 64;i++){
        a = this->polyarray[i];
        res->polyarray[i] = -this->polyarray[128-i];
        res->polyarray[128-i] = -a;
    }

}
template<typename T>
void Poly<T>::copy(Poly * target){
    this->nttflag = target->nttflag;
    for(int i = 0 ; i < 128; i++){
        this->polyarray[i] = target->polyarray[i];
    }
}
template <typename T>
T Poly<T>::mod(T a ,int flag){
    // flag = 0: [-q/2,q/2] , 1:[0,q]
    T res = a % q;
    if(flag ==0){
        if(res > q/2){
            return res - q;
        }
        if(res < -q/2){
            return res + q;
        }
    }
    else{
        if(res < 0){
            return res + q;
        }
    }
    return res;
}
template <typename T>
T Poly<T>::mulmod(T a, T b)
{
    T res = 0;
    int flag;
    a = mod(a);
    b = mod(b);
    while(b != 0){
        if((b&1 != 0)){
            res = mod(res + a);
        }
        b >>= 1;
        a = mod(a << 1);
    }
    return res;
}
template <typename T>
void Poly<T>::output(int flag)
{
    if(nttflag == 1){
        std::cout << "nttarray: ";
    }
    else{
        std::cout << "polyarray: ";
    }
    for(int i = 0 ; i < 128;i++){
        std::cout << mod(polyarray[i] , flag) << " ";
    }
    std::cout << std::endl;
}

template <typename T>
int64_t Poly<T>::norm(){
    this->to_poly();
    int64_t res = 0, temp;
    for(int i = 0 ; i < 128;i++){
        temp = mod(polyarray[i] , 0);
        res += temp*temp;
    }
    return res;
}

template class Poly<int16_t>;
template class Poly<int64_t>;


int64_t NTTTABLEQ[128] = {1,-19713449676,-12684577924,-23224220986,-22575945377,-32427681721, 7201387,3120114039,
    13416529352,-32119515348,-27584096198,-33356371947,10297361342,-25435822379,-28928468348,-1213100144,
    -25367629841,892035260,-31598040742,24253212821,11912019664,26423810131,-1883931092,-4818611860,
    2479314337,-29811804473,-27267549940,7169898790,-18420364768,-15547629305,15809916002,13863862447,
    2937993509,-5385229747,-33518795608,28827130605,6909479775,-11859717328,3329469899,21531097251,
    4130111020,-29495453243,26052388949,25183604338,-17974785716,-23139277071,-14203244871,-6330672295,
    4759051306,-22959701434,-34006657221,-27138458928,18185935589,-9294923434,-3273073298,26563800281,
    24350991323,3522840960,784954030,16109869091,-26382811400,31682606952,29243827564,-31363260053,
    -9694590617,-20942563199,-1749249764,-9796659360,-1616107111,10035157443,14696758156,23433501642,
    17934269046,-3994019765,-23487210416,-23123015975,-340392228,-26245707274,-32580107802,-17066224306,
    -20679227766,21088720345,-1504163901,31627925662,-25793757978,-22653835460,21336314817,3353503744,
    -30198214964,-12519193212,14251865123,15293499961,13101246360,-17626455909,-20172529677,-10226476308,
    -18405005580,-11508236726,-21933923436,7948321023,-16221629165,-17114748789,8019688474,705654232,
    24236819565,-5625856111,19348262037,-32711297059,12475839286,33434844938,-6045265618,8086816399,
    -14563246824,-30055451481,-14609432832,-30791895383,29592085416,21390685797,22479282053,31664758284,
    20987596018,30727804836,-22593640817,-32066643745,26915213124,23189600899,-11849453608,8712564243,
        };

Polyq::Polyq(int64_t *array, int nttflag):Poly(array , nttflag){
    q = Q;
    q2inv = 68182597951;
    ntttable = NTTTABLEQ;
}

Polyq::Polyq(int nttflag):Poly(nttflag){
    q = Q;
    q2inv = 68182597951;
    ntttable = NTTTABLEQ;
}


Polyq::Polyq(unsigned char *seed, int seedlen, int eta, int nttflag):Poly(nttflag){
    // eta = 0: uniform ; eta != 0 : cbd with param eta
    q = Q;
    q2inv = 68182597951;
    ntttable = NTTTABLEQ;
    keccak_state state;
    shake256_init(&state);
    shake256_absorb(&state , seed , seedlen);
    
    if(eta == 0){
        uniformq(&state , this->polyarray);
    }
    else{
        cbd<int64_t>(&state , this->polyarray , eta);
    }
}

Polyq::Polyq(keccak_state *state, int eta, int nttflag):Poly(nttflag){
    q = Q;
    q2inv = 68182597951;
    ntttable = NTTTABLEQ;
    if(eta == 0){
        uniformq(state , this->polyarray);
    }
    else{
        cbd<int64_t>(state , this->polyarray , eta);
    }
}

Polyq::Polyq(double sigma, int nttflag):Poly(nttflag){
    q = Q;
    q2inv = 68182597951;
    ntttable = NTTTABLEQ;
    Gaussian<int64_t>(this->polyarray , sigma);
}

void Polyq::from_char(unsigned char *polychar, int nttflag, int eta)
{
    this->nttflag = nttflag;
    ino64_t c[8];
    if(eta == 0){
        for(int i = 0;i<64;i++){
            polyarray[i*2] = 0;
            polyarray[i*2+1] = 0;
            for(int j = 0 ; j < 4; j++){
                polyarray[i*2] <<= 8;
                polyarray[i*2] |= polychar[i*9+j];
            }
            polyarray[i*2] <<= 4;
            polyarray[i*2] |= polychar[i*9+4]>>4;
            polyarray[i*2+1] = polychar[i*9+4]&0xf;
            for(int j = 5 ; j < 9; j++){
                polyarray[i*2+1] <<= 8;
                polyarray[i*2+1] |= polychar[i*9+j];
            }
        }
    }
    else{
        assert(eta == 7);
        for(int i = 0; i< 64; i++){
            polyarray[i*2] = (polychar[i]>>4)-eta;
            polyarray[i*2+1] = (polychar[i]&0xf) - eta;
        }
    }
}

void Polyq::to_char(unsigned char *res , int eta){
    int64_t a,b;
    int64_t c[8];
    if(eta == 0){
        for(int i = 0 ;i<64;i++){
            a = mod(polyarray[i*2]);
            b = mod(polyarray[i*2+1]);
        
            for(int j = 0 ; j < 4; j++){
                res[i*9+j] = (a>>(8*(3-j) + 4))&0xff;
            }
            res[i*9+4] = ((a&0xf)<<4) | (b>>32);
            for(int j = 5;j < 9;j++){
                res[i*9+j] = (b>>(8*(8-j)))&0xff;
            }
        }
    }
    else{
        assert(eta == 7);
        this->to_poly();

        for(int i = 0 ; i < 64 ;i++){
            res[i] = ((polyarray[2*i]+eta)<<4)|(polyarray[2*i+1]+eta);
        }
    }
    return;
}
int16_t ntt3329[] = {
       1,   -1600,    -749,     -40,    -687,     630,   -1432,     848,
    1062,   -1410,     193,     797,    -543,     -69,     569,   -1583,
     296,    -882,    1339,    1476,    -283,      56,   -1089,    1333,
    1426,   -1235,     535,    -447,    -936,    -450,   -1355,     821,
     289,     331,     -76,   -1573,    1197,   -1025,   -1052,   -1274,
     650,   -1352,    -816,     632,    -464,      33,    1320,   -1414,
   -1010,    1435,     807,     452,    1438,    -461,    1534,    -927,
    -682,    -712,    1481,     648,    -855,    -219,    1227,     910,
      17,    -568,     583,    -680,    1637,     723,   -1041,    1100,
    1409,    -667,     -48,     233,     756,   -1173,    -314,    -279,
   -1626,    1651,    -540,   -1540,   -1482,     952,    1461,    -642,
     939,   -1021,    -892,    -941,     733,    -992,     268,     641,
    1584,   -1031,   -1292,    -109,     375,    -780,   -1239,    1645,
    1063,     319,    -556,     757,   -1230,     561,    -863,    -735,
    -525,    1092,     403,    1026,    1143,   -1179,    -554,     886,
   -1607,    1212,   -1455,    1029,   -1219,    -394,     885,   -1175,
};

Polyp::Polyp(int16_t *array, int nttflag):Poly(array , nttflag){
    q = 3329;
    q2inv = 3303;
    ntttable = ntt3329;
}

Polyp::Polyp(int nttflag):Poly(nttflag){
    q = 3329;
    q2inv = 3303;
    ntttable = ntt3329;
}


Polyp::Polyp(unsigned char *seed, int seedlen, int eta, int nttflag):Poly(nttflag){
    // eta = 0: uniform ; eta != 0 : cbd with param eta
    q = 3329;
    q2inv = 3303;
    ntttable = ntt3329;
    keccak_state state;
    shake256_init(&state);
    shake256_absorb(&state , seed , seedlen);
    
    if(eta == 0){
        uniformp(&state , this->polyarray);
    }
    else{
        cbd<int16_t>(&state , this->polyarray , eta);
    }
}

Polyp::Polyp(keccak_state *state, int eta, int nttflag):Poly(nttflag){
    q = 3329;
    q2inv = 3303;
    ntttable = ntt3329;
    if(eta == 0){
        uniformp(state , this->polyarray);
    }
    else{
        cbd<int16_t>(state , this->polyarray , eta);
    }
}
Polyp::Polyp(double sigma, int nttflag):Poly(nttflag){
    q = 3329;
    q2inv = 3303;
    ntttable = ntt3329;
    Gaussian(this->polyarray,sigma);
}
void Polyp::from_char(unsigned char *polychar, int nttflag, int eta)
{
    this->nttflag = nttflag;
    int16_t c[8];
    if(eta == 0){
        for(int i = 0;i<64;i++){
            polyarray[2*i] = polychar[i*3]<< 4;
            polyarray[2*i] |= polychar[i*3+1]>>4;
            polyarray[2*i+1] = (polychar[i*3+1]&0xf)<<8;
            polyarray[2*i+1] |= polychar[i*3+2];
            polyarray[2*i] = mod(polyarray[2*i]);
            polyarray[2*i+1] = mod(polyarray[2*i+1]);
        }
    }
    else{
        assert(eta == 7);
        for(int i = 0; i< 64; i++){
            polyarray[i*2] = (polychar[i]>>4)-eta;
            polyarray[i*2+1] = (polychar[i]&0xf) - eta;
        }
    }
}

void Polyp::to_char(unsigned char *res , int eta){
    int16_t a,b;
    int16_t c[8];
    if(eta == 0){
        for(int i = 0 ;i<64;i++){
            a = mod(polyarray[i*2]);
            b = mod(polyarray[i*2+1]);
        
            res[i*3] = a>>4;
            res[i*3+1] = ((a&0xf)<<4) | (b>>8);
            res[i*3+2] = b&0xff;
        }
    }
    else{
        assert(eta == 7);
        this->to_poly();

        for(int i = 0 ; i < 64 ;i++){
            res[i] = ((polyarray[2*i]+eta)<<4)|(polyarray[2*i+1]+eta);
        }
    }
    return;
}

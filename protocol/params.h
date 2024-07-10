#define polymatp polymat<Polyp>
#define polymatq polymat<Polyq>
#define polyvecp polyvec<Polyp>
#define polyvecq polyvec<Polyq>


#define N1 4
#define M2 13
#define M 2
#define M1 (3*N+M+K1+3)
#define TAU 4
#define BETAE 81536
#define BETAR 1536
#define BETAE2 6272
#define T 1.64
#define Y1 7376.40
#define Y2 4283.18
#define Y3 2868.60
#define Y4 1128103.54
#define B1 817683.17
#define B2 247091.6
#define B3 (T*842567.93)
#define B4 289929193.97

#define PKLEN (VECLEN+512)
#define SKLEN (VECETALEN*4+256)
#define PI1LEN ((13+2*TAU+N1*2+N*5+M*2+K1*2)*POLYLEN+16*256)
#define PI2LEN ((16+TAU*2+2*N1+6*N)*POLYLEN+8*256)
#define UPOLYLEN (5*32)
#define VPOLYLEN 64
#define UVECLEN (UPOLYLEN*4)
#define polymatp polymat<Polyp>
#define polymatq polymat<Polyq>
#define polyvecp polyvec<Polyp>
#define polyvecq polyvec<Polyq>


#define N1 4
#define N_BATCH 10
#define M2 13
#define M 2
#define M1 (3*N+M+K1+3)
#define M1_2 (2*N+6)
#define M1_BATCH (4*N+6+N_BATCH)
#define TAU 4
#define BETAE 81536
#define BETAR 1536
#define BETAE2 6272
#define NBETAE2 (N_BATCH*BETAE2)
#define T 1.70
#define Y1 7376.40
#define Y2 4283.18
#define Y3 2868.60*10
#define Y4 1128103
#define B1 (817683.17)
#define B2 247091.6
#define B3 ((T*842567.93)*1.5)
#define B4 (289929193.97 * 128)

#define B11 (817683.17*1.5)
#define B22 247091.6
#define B33 ((T*842567.93)*1.5)

#define PKLEN (2*VECLEN+32)
#define SKLEN (VECETALEN*2+32)
#define PI1LEN ((12+2*TAU+N1*2+N*3+M+K1+M2)*POLYLEN+16*256)
#define PI2LEN ((16+TAU*2+2*N1+6*N)*POLYLEN+8*256)
#define PIBLEN (POLYLEN*(4+TAU+2*N1 + M1_BATCH+M2) + 32+VECETALEN + 256*8 + TAU*POLYLEN + POLYLEN*2)
#define UPOLYLEN (5*32)
#define VPOLYLEN 64
#define UVECLEN (UPOLYLEN*4)
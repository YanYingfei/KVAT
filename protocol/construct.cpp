#include "construct.h"
#include <iostream>

void constructFy(polymatq *Fymat ,polymatq *F1 ,polymatq *Fg ,polymatq *F3 ,polymatq *F4){
    for(int i = 0 ; i < N ; i++){
        Fymat->vecarray[i] = F1->vecarray[i];
    }
    for(int i = 0 ; i < TAU ; i++){
        Fymat->vecarray[N+i] = Fg->vecarray[i];
    }
    for(int i = 0 ; i < 2;i++){
        Fymat->vecarray[N+TAU+i] = F3->vecarray[i];
    }
    for(int i = 0 ; i < 2;i++){
        Fymat->vecarray[N+TAU+2+i] = F4->vecarray[i];
    }
}
void constructy(polyvecq *y ,polymatq *Fymat ,polyvecq *y1 ,polyvecq*y2){
    polyvecq temp(y1->k);
    y1->sigma(&temp);
    for(int i = 0 ; i < y1->k ; i++){
        y->polyarray[i]->copy(y1->polyarray[i]);
        y->polyarray[i+y1->k]->copy(temp.polyarray[i]);
    }
    polyvecq temp1(Fymat->k);
    Fymat->right_mul(&temp1 , y2);
    polyvecq temp2(Fymat->k);
    temp1.sigma(&temp2);
    for(int i = 0 ; i < temp1.k ; i++){
        temp1.polyarray[i]->mul_num(y->polyarray[i+y1->k*2] , -1);
        temp2.polyarray[i]->mul_num(y->polyarray[i+y1->k*2+temp1.k] , -1);
    }

}
void constructD2(Polyq *D2list[] ,polyvecq *mu ,polymatq *Gamma){
    
    Polyq temp(1);
    for(int i = 0 ; i < 4;i++){
        D2list[i] = new Polyq(1);
        for(int j = 0 ; j < TAU ; j++){
            mu->polyarray[j]->mul(&temp , Gamma->vecarray[j]->polyarray[N+i]);
            D2list[i]->add(D2list[i] , &temp);
        }
    }
    Polyq *temp2;
    temp2 = D2list[0];
    D2list[0] = D2list[1];
    D2list[1] = D2list[2];
    D2list[2] = temp2;
    D2list[2]->mul_num(D2list[2] , -1);

}
void constructd1(polyvecq *d1vec ,polyvecq* mu ,polymatq *Gamma ,polymatq *A){
    
    polyvecq temp1(N , 1) , temp2(N , 1);
    polyvecq temp3(N , 1) , temp4(N , 1);

    for(int j = 0 ; j < TAU ; j++){
        temp2.reset(1);
        for(int i = 0 ; i < N ; i++){
            A->vecarray[i]->mul_poly(&temp1 , Gamma->vecarray[j]->polyarray[i]);
            temp2.add(&temp2 , &temp1);
        }
        temp2.mul_poly(&temp3 , mu->polyarray[j]);
        temp4.add(&temp4 , &temp3);
    }
    for(int i = 0 ; i < N ; i++){
        d1vec->polyarray[i]->copy(temp4.polyarray[i]);
    }

    
    temp2.reset(1);

    for(int j = 0 ; j < TAU ; j++){
        for(int i = 0 ; i < N ; i++){
            temp1.polyarray[i]->copy(Gamma->vecarray[j]->polyarray[i]);
        }
        temp1.mul_poly(&temp1 , mu->polyarray[j]);
        temp2.add(&temp2 , &temp1);
    }
    for(int i = 0 ; i < N ; i++){
        d1vec->polyarray[i+N+1]->copy(temp2.polyarray[i]);
        d1vec->polyarray[i+2*(2*N+M+K1+3)]->copy(temp2.polyarray[i]);
    }

    Polyq temp5(1), temp6(1);
    Polyq temp7(0);
    for(int i = 0 ; i < 128;i++){
        temp7.polyarray[i] = 1;
    }
    temp7.sigma(&temp7);
    temp6.reset(1);
    for(int j = 0 ; j < TAU ; j++){
        Gamma->vecarray[j]->polyarray[N]->mul(&temp5 , mu->polyarray[j]);
        
        temp6.add(&temp6 , &temp5);
    }
    temp6.mul(&temp6 , &temp7);
    d1vec->polyarray[2*N+2]->copy(&temp6);
    d1vec->polyarray[2*N+3]->copy(&temp6);

    for(int j = 0 ; j < TAU ; j++){
        d1vec->polyarray[j+2*(2*N+M+K1+3)+N]->copy(mu->polyarray[j]);
    }

}

void constructm(polyvecq *m ,polyvecq *m1 ,polyvecq *m2 ,polyvecq *g ,polyvecq *y3 ,polyvecq *y4){
    int mlen = 2*N+M+K1+3;
    int m2len = N+TAU+4;
    for(int i = 0 ; i < mlen;i++){
        m->polyarray[i]->copy(m1->polyarray[i]);
        m1->polyarray[i]->sigma(m->polyarray[i+mlen]);
    }
    for(int i = 0 ; i < N; i++){
        m->polyarray[2*mlen+i]->copy(m2->polyarray[i]);
        m2->polyarray[i]->sigma(m->polyarray[i+mlen*2+m2len]);
    }
    for(int i = 0 ; i < TAU; i++){
        m->polyarray[2*mlen+N+i]->copy(g->polyarray[i]);
        g->polyarray[i]->sigma(m->polyarray[i+mlen*2+m2len+N]);
    }
    for(int i = 0 ; i < 2; i++){
        m->polyarray[2*mlen+N+TAU+i]->copy(y3->polyarray[i]);
        y3->polyarray[i]->sigma(m->polyarray[i+mlen*2+m2len+N+TAU]);
    }
    for(int i = 0 ; i < 2; i++){
        m->polyarray[2*mlen+N+TAU+2+i]->copy(y4->polyarray[i]);
        y4->polyarray[i]->sigma(m->polyarray[i+mlen*2+m2len+N+TAU+2]);
    }
}

void constructz(polyvecq *z ,polymatq*Fymat, polyvecq *tF, polyvecq *tg, polyvecq *t3, polyvecq *t4,Polyq*cpoly, polyvecq *z1, polyvecq *z2){
    polyvecq ty(TAU+N+4);
    for(int i = 0 ; i < N; i++){
        ty.polyarray[i]->copy(tF->polyarray[i]);
    }
    for(int i = 0 ; i < TAU; i++){
        ty.polyarray[i+N]->copy(tg->polyarray[i]);
    }
    for(int i = 0 ; i < 2; i++){
        ty.polyarray[i+N+TAU]->copy(t3->polyarray[i]);
    }
    for(int i = 0 ; i < 2; i++){
        ty.polyarray[i+N+TAU+2]->copy(t4->polyarray[i]);
    }
    polyvecq temp(z1->k);
    z1->sigma(&temp);
    for(int i = 0 ; i < z1->k ; i++){
        z->polyarray[i]->copy(z1->polyarray[i]);
        z->polyarray[i+z1->k]->copy(temp.polyarray[i]);
    }
    polyvecq temp1(Fymat->k);
    Fymat->right_mul(&temp1 , z2);
    polyvecq temp2(Fymat->k);
    ty.mul_poly(&ty , cpoly);
    ty.sub(&temp1 , &temp1);
    temp1.sigma(&temp2);
    for(int i = 0 ; i < temp1.k ; i++){
        z->polyarray[i+z1->k*2]->copy(temp1.polyarray[i]);
        z->polyarray[i+z1->k*2+temp1.k]->copy(temp2.polyarray[i]);
    }
}

void constructd0(Polyq *d, polyvecq *mu, polymatq *Gamma, polyvecq *c , polyvecq *h){
    Polyq temp1(1) , temp2(1);
    Polyq temp3(1) , temp4(1);
    d->to_ntt();
    for(int j = 0 ; j < TAU ; j++){
        temp2.reset(1);
        for(int i = 0 ; i < N ; i++){
            c->polyarray[i]->mul(&temp1 , Gamma->vecarray[j]->polyarray[i]);
            temp2.add(&temp2 , &temp1);
        }
        temp2.mul(&temp3 , mu->polyarray[j]);
        d->add(d , &temp3);
    }
    d->to_poly();

    for(int j = 0 ; j < TAU; j++){
        temp1.reset(0);
        Gamma->vecarray[j]->polyarray[N+1]->add(&temp1 , Gamma->vecarray[j]->polyarray[N+2]);
        temp1.mul(&temp2 , mu->polyarray[j]);
        temp2.mul_num(&temp2 , BETAE);
        d->add(d , &temp2);
    }

    for(int j = 0 ; j < TAU; j++){
        Gamma->vecarray[j]->polyarray[N+3]->mul(&temp2 , mu->polyarray[j]);
        temp2.mul_num(&temp2 , BETAR);
        d->add(d , &temp2);
    }

    h->mul(&temp1 , mu);
    temp1.to_poly();
    d->add(d , &temp1);
    d->mul_num(d , -1);

    
    
}

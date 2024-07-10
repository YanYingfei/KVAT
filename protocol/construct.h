#ifndef CONSTRUCT_H
#define CONSTRUCT_H

#include "../poly/polymat.h"
#include "../poly/poly.h"
#include "params.h"



void constructFy(polymatq *Fymat ,polymatq *F1 ,polymatq *Fg ,polymatq *F3 ,polymatq *F4);
void constructy(polyvecq *y ,polymatq *Fymat ,polyvecq *y1 ,polyvecq*y2);
void constructD2(Polyq *D2list[] ,polyvecq *mu ,polymatq *Gamma);
void constructd1(polyvecq *d1vec ,polyvecq* mu ,polymatq *Gamma ,polymatq *A);
void constructm(polyvecq *m ,polyvecq *m1 ,polyvecq *m2 ,polyvecq *g ,polyvecq *y3 ,polyvecq *y4);
void constructz(polyvecq *z ,polymatq*Fymat, polyvecq *tE, polyvecq *tg, polyvecq *t3, polyvecq *t4,Polyq*cpoly, polyvecq *z1, polyvecq *z2);
void constructd(Polyq *d ,polyvecq* mu ,polymatq* Gamma ,polyvecq* c);
#endif
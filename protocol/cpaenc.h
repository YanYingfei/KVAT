#ifndef KYBER_H
#define KYBER_H

void keygen(unsigned char* pke);
void enc(unsigned char *cipher,unsigned char* uv , unsigned char* rseed, unsigned char *pke, unsigned char*m , int mlen);

#endif
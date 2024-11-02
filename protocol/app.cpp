#include "main_protocol.h"
#include <iostream>



int main(){
    unsigned char pk[PKLEN] , sk[SKLEN];
    unsigned char query[VECLEN],pi1[PI1LEN] , crsseed[32],H1input[32+MD/8] , xchar[VECETALEN];
    unsigned char resp[VECLEN] , pi2[PI2LEN] , crsseed2[32];
    unsigned char msg[VECLEN+32+MD/8];
    int bit;
    
    keygen(pk , sk);

    AT_CQ(pi1 ,crsseed, query , H1input , xchar , pk);
    
    AT_IT(pi2 , crsseed2 , resp , pk , sk , query , crsseed , pi1);
    AT_CF(msg , pk , resp , H1input , xchar , query , pi2 , crsseed2);
    bit = AT_Rb(pk , sk , msg);
    std::cout << bit << std::endl;
}
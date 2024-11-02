#include "main_protocol.h"
#include <ctime>
#include <iostream>
double max(double list[] , int len){
    double res = list[0];
    for(int i = 0 ; i < len ; i++){
        if(list[i] > res){
            res = list[i];
        }
    }
    return res;
}
double min(double list[] , int len){
    double res = list[0];
    for(int i = 0 ; i < len ; i++){
        if(list[i] < res){
            res = list[i];
        }
    }
    return res;
}
double average(double list[] , int len){
    double res = 0;
    for(int i = 0 ; i < len ; i++){
        res += list[i];
    }
    res /= len;
    return res;
}
int main(){
    unsigned char pk[PKLEN] , sk[SKLEN];
    unsigned char query[VECLEN],pi1[PI1LEN] , crsseed[32],H1input[32+MD/8] , xchar[VECETALEN];
    unsigned char resp[VECLEN] , pi2[PI2LEN] , crsseed2[32];
    unsigned char msg[VECLEN+32+MD/8];
    int n = 5;
    double kgtime[n] , CQtime[n] , ITtime[n] , CFtime[n] ,\
    Rbtime[n], pi1prvtime[n] , pi1vrftime[n],pi2prvtime[n] , pi2vrftime[n];
    time_t t1 , t2;
    for(int i = 0 ; i < n ; i++){
        std::cout << " \n time: "<< i << std::endl;
        t1 = clock();
        keygen(pk , sk);
        t2 = clock();
        kgtime[i] = (double)(t2-t1)/CLOCKS_PER_SEC;

        t1 = clock();
        AT_CQ(pi1 ,crsseed, query , H1input , xchar , pk , pi1prvtime+i);
        t2 = clock();
        CQtime[i] = (double)(t2-t1)/CLOCKS_PER_SEC;

        t1 = clock();
        AT_IT(pi2 , crsseed2 , resp , pk , sk , query , crsseed , pi1 ,pi1vrftime+i ,pi2prvtime+i);
        t2 = clock();
        ITtime[i] = (double)(t2-t1)/CLOCKS_PER_SEC;

        t1 = clock();
        AT_CF(msg , pk , resp , H1input , xchar , query , pi2 , crsseed2 , pi2vrftime+i);
        t2 = clock();
        CFtime[i] = (double)(t2-t1)/CLOCKS_PER_SEC;

        t1 = clock();
        AT_Rb(pk , sk , msg);
        t2 = clock();
        Rbtime[i] = (double)(t2-t1)/CLOCKS_PER_SEC;
    }
    std::cout << max(kgtime , n) <<" "<< min(kgtime , n) << " " << average(kgtime , n) << std::endl;
    std::cout << max(CQtime , n) <<" "<< min(CQtime , n) << " " << average(CQtime , n) << std::endl;
    std::cout << max(ITtime , n) <<" "<< min(ITtime , n) << " " << average(ITtime , n) << std::endl;
    std::cout << max(CFtime , n) <<" "<< min(CFtime , n) << " " << average(CFtime , n) << std::endl;
    std::cout << max(Rbtime , n) <<" "<< min(Rbtime , n) << " " << average(Rbtime , n) << std::endl;
    std::cout << max(pi1prvtime , n) <<" "<< min(pi1prvtime , n) << " " << average(pi1prvtime , n) << std::endl;
    std::cout << max(pi1vrftime , n) <<" "<< min(pi1vrftime , n) << " " << average(pi1vrftime , n) << std::endl;
    std::cout << max(pi2prvtime , n) <<" "<< min(pi2prvtime , n) << " " << average(pi2prvtime , n) << std::endl;
    std::cout << max(pi2vrftime , n) <<" "<< min(pi2vrftime , n) << " " << average(pi2vrftime , n) << std::endl;
    return 0;
}
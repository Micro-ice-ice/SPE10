#pragma once
#include <cmath>
#include "vars.hpp"

inline double K_si(double pb, double pi, double si, bool w){
    if (pb > pi){
       if (w){
           return 1;
       }
       else return 0;
    }
    else{
        if (!w){
            return 1 - si;
        }
        else return si;
    }
}

inline double WI(double kx, double ky){
    return (2 * PI * HZ * sqrt(kx * ky)) /
            log(0.28 * sqrt(pow(kx, 2)
            * sqrt(ky / kx) + pow(ky, 2) * sqrt(kx / ky))
            / (RW * pow(ky / kx, 0.25) + pow(kx / ky, 0.25) +
            S));
}

inline double T_ij(double ki, double kj, double h){
    return 2 * ki * kj / (h * (ki + kj));
}

// inline double T_boundary(double ki, double h, double u_init = 1) {
//     return 2 * ki * u_init / pow(h, 2);
// }

inline double IfFuncRp(double si, double sj, double pi, double pj){
    if (pi > pj){
        return si;
    } else return sj;
}

inline double IfFuncRs(double pi, double pj){
    if (pi > pj){
        return 1;
    } else return 0;
}

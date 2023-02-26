#pragma once
#include <cmath>
#include "vars.hpp"
// только для dR/dpi и dR/dpj
inline double K_p(double pb, double pi, double si){
    if (pb > pi){
        return 0;
    }
    else return si;
}

// для dR0/dsi и dRw/dsi (для элементов вне диагонали не используется)
inline double K_s(double pb, double pi, double si, bool w){
    if (!w){
        if (pb > pi){
            return 0;
        }
        else return 1;
    }
    else {
        if (pb > pi) {
            return 1;
        } else return si;
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

inline double IfFuncRpi(double si, double sj_top, double sj_down, double sj_left, double sj_right,
                        double pi, double pj_top, double pj_down, double pj_left, double pj_right){
    if (pi < pj_top) {
        return sj_top;
    }
    if (pi < pj_down) {
        return sj_down;
    }
    if (pi < pj_left) {
        return sj_left;
    }
    if (pi < pj_right) {
        return sj_right;
    }
    else return si;
}

inline double IfFuncRs(double pi, double pj, bool diag){
    if (diag) {
        if (pi > pj) {
            return 1;
        } else return 0;
    }
    else{
        if (pi > pj) {
            return 0;
        } else return 1;
    }
}


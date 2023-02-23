#pragma once
#include <cmath>
#include "vars.hpp"

double k_si(double pb, double pi, double Si, bool w){
    if (pb > pi){
       if (w){
           return 1;
       }
       else return 0;
    }
    else{
        if (!w){
            return 1 - Si;
        }
        else return Si;
    }
} //

double WI(double kx, double ky){
    return (2 * PI * HZ * sqrt(kx * ky)) /
            log(0.28 * sqrt(pow(kx, 2)
            * sqrt(ky / kx) + pow(ky, 2) * sqrt(kx / ky))
            / (RW * pow(ky / kx, 0.25) + pow(kx / ky, 0.25) +
            S));
}


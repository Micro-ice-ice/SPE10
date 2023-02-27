#pragma once
#include <cmath>
#include <vector>
#include "vars.hpp"
#include "solver/pcg.cpp"

using namespace std;

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

inline vector<double> Solve(vector<int> jacobi_row, vector<int> jacobi_column, vector<double> jacobi_value, vector<double> r){
    //write J.mtx, R.txt
    std::ofstream output("J.mtx"),  vec("r.txt");
    output << "%%MatrixMarket matrix coordinate real general" << std::endl;
    output << 2 * NX * NY << " " << 2 * NX * NY  << " " << jacobi_value.size() << std::endl;
    vec << 2 * NX * NY << std::endl;

    for (size_t i = 0; i < jacobi_row.size() - 1; ++i)
    {
        for (int j = jacobi_row[i]; j < jacobi_row[i + 1]; ++j)
            output << i + 1 << " " << jacobi_column[j] + 1 << " " << jacobi_value[j] << std::endl;
        vec << r[i] << std::endl;

    }
    output.close();

    system("pcg.exe");

    string filename = "../solution_pcg.txt";   // Name of the file
    ifstream newfile (filename);
    vector<double>delta_R;
    if (newfile.is_open()){ //checking whether the file is open
        //cout << "File is open";
        double b;
        while (newfile >> b){
            delta_R.push_back(b);
        }
        newfile.close(); //close the file object.
    } else {
        cout << "File isn't open";
    }

    return delta_R;
}
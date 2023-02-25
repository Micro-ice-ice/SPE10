#include <iostream>
#include <fstream>
#include <string>
#include "vars.hpp"
#include "f.hpp"
#include "cell.hpp"
#include <vector>

using namespace std;

int main(int argc, char **argv) {

    Cell::SetNx(NX);
    Cell::SetNy(NY);
    Cell::SetCellCount(NX * NY);

    vector<Cell> cells;

    double *grid = new double[NX * NY * 3];
    string filename = "../por_perm_case2a/data.txt";   // Name of the file
    ifstream newfile (filename);
    if (newfile.is_open()){ //checking whether the file is open
        //cout << "File is open";
        int i = 0;
        double kx, ky, phi;
        while (newfile >> kx >> ky >> phi) {

            cells.push_back(Cell(kx, ky, phi));
        }
        newfile.close(); //close the file object.
    } else {
        cout << "File isn't open";
    }

    return 0;

}

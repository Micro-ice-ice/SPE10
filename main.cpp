#include <iostream>
#include <fstream>
#include <string>
#include "vars.hpp"
#include "f.hpp"

using namespace std;

int main(int argc, char **argv) {

    double *grid = new double[NX * NY * 3];
    string filename = "../por_perm_case2a/data.txt";   // Name of the file
    ifstream newfile (filename);
    if (newfile.is_open()){ //checking whether the file is open
        //cout << "File is open";
        int i = 0;
        double value;
        while (newfile >> value) {

            grid[i] = value;
            ++i;
        }
        newfile.close(); //close the file object.
    } else {
        cout << "File isn't open";
    }

    return 0;

}

#include <iostream>
#include <fstream>
#include <string>
#include "vars.hpp"
#include "f.hpp"
#include "cell.hpp"
#include <vector>

using namespace std;

int main(int argc, char **argv) {

    //init grid parameters

    Cell::SetNx(NX);
    Cell::SetNy(NY);
    Cell::SetCellCount(NX * NY);

    vector<Cell> cells; //grid

    vector<double> jacobi_value; //COO jacobi mattrix format
    vector<int> jacobi_row;
    vector<int> jacobi_column;

    //open data file

    string filename = "../por_perm_case2a/data.txt";   // Name of the file
    ifstream newfile (filename);
    if (newfile.is_open()){ //checking whether the file is open
        //cout << "File is open";
        double kx, ky, phi;
        while (newfile >> kx >> ky >> phi) {

            cells.push_back(Cell(kx, ky, phi));
        }
        newfile.close(); //close the file object.
    } else {
        cout << "File isn't open";
    }

    //create jacobi

    for(int i = 0; i < cells.size(); ++i){

        Cell &cell = cells[i];

        //first block 
        {
            double sum = 0;

            //top 

            if (cell.Top == true){

                Cell &topCell = cells[cell.GetTopIndex()];

                double topValue = 
                T_ij(cell.GetKy(), topCell.GetKy(), HY) * 
                IfFuncRp(cell.GetS(), topCell.GetS(), cell.GetP(), topCell.GetP()) + 
                K_si(0, cell.GetP(), cell.GetS(), false) * 
                WI(cell.GetKx(), cell.GetKy());

                sum = sum + topValue;
            }

            

            //bottom

            //left

            //right

        }
    }

    return 0;

}

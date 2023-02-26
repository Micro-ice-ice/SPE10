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

        //first block oil dR/dP
        {
            double sum = 0;

            //top 

            if (cell.Top == true){

                Cell &topCell = cells[cell.GetTopIndex()];

                double topValue = 
                T_ij(cell.GetKy(), topCell.GetKy(), HY) * 
                IfFuncRp(cell.GetS(), topCell.GetS(), cell.GetP(), topCell.GetP()) + 
                K_si(PB, cell.GetP(), cell.GetS(), false) * 
                WI(cell.GetKx(), cell.GetKy());

                jacobi_value.push_back(-topValue);
                jacobi_row.push_back(i);
                jacobi_column.push_back(cell.GetTopIndex());

                sum = sum + topValue;
            }

            //bottom

            if (cell.Bottom == true){

                Cell &bottomCell = cells[cell.GetBottomIndex()];

                double bottomValue = 
                T_ij(cell.GetKy(), bottomCell.GetKy(), HY) * 
                IfFuncRp(cell.GetS(), bottomCell.GetS(), cell.GetP(), bottomCell.GetP()) + 
                K_si(PB, cell.GetP(), cell.GetS(), false) * 
                WI(cell.GetKx(), cell.GetKy());

                jacobi_value.push_back(-bottomValue);
                jacobi_row.push_back(i);
                jacobi_column.push_back(cell.GetBottomIndex());

                sum = sum + bottomValue;
            }

            //left

            if (cell.Left == true){

                Cell &leftCell = cells[cell.GetLeftIndex()];

                double leftValue = 
                T_ij(cell.GetKx(), leftCell.GetKx(), HX) * 
                IfFuncRp(cell.GetS(), leftCell.GetS(), cell.GetP(), leftCell.GetP()) + 
                K_si(PB, cell.GetP(), cell.GetS(), false) * 
                WI(cell.GetKx(), cell.GetKy());

                jacobi_value.push_back(-leftValue);
                jacobi_row.push_back(i);
                jacobi_column.push_back(cell.GetLeftIndex());

                sum = sum + leftValue;
            }

            //right

             if (cell.Right == true){

                Cell &rightCell = cells[cell.GetRightIndex()];

                double rightValue = 
                T_ij(cell.GetKx(), rightCell.GetKx(), HX) * 
                IfFuncRp(cell.GetS(), rightCell.GetS(), cell.GetP(), rightCell.GetP()) + 
                K_si(PB, cell.GetP(), cell.GetS(), false) * 
                WI(cell.GetKx(), cell.GetKy());

                jacobi_value.push_back(-rightValue);
                jacobi_row.push_back(i);
                jacobi_column.push_back(cell.GetRightIndex());

                sum = sum + rightValue;
            }

            //cell 

            jacobi_value.push_back(sum);
            jacobi_row.push_back(i);
            jacobi_column.push_back(i);

        }

        //firth block water dR/dP
        {
            int shift = NX * NY;
            double sum = 0;

            //top 

            if (cell.Top == true){

                Cell &topCell = cells[cell.GetTopIndex()];

                double topValue = 
                T_ij(cell.GetKy(), topCell.GetKy(), HY) * 
                IfFuncRp(1 - cell.GetS(), 1 - topCell.GetS(), cell.GetP(), topCell.GetP()) + 
                K_si(PB, cell.GetP(), 1 - cell.GetS(), true) * 
                WI(cell.GetKx(), cell.GetKy());

                jacobi_value.push_back(-topValue);
                jacobi_row.push_back(i + shift);
                jacobi_column.push_back(cell.GetTopIndex() + shift);

                sum = sum + topValue;
            }

            //bottom

            if (cell.Bottom == true){

                Cell &bottomCell = cells[cell.GetBottomIndex()];

                double bottomValue = 
                T_ij(cell.GetKy(), bottomCell.GetKy(), HY) * 
                IfFuncRp(1 - cell.GetS(), 1 - bottomCell.GetS(), cell.GetP(), bottomCell.GetP()) + 
                K_si(PB, cell.GetP(), 1 - cell.GetS(), true) * 
                WI(cell.GetKx(), cell.GetKy());

                jacobi_value.push_back(-bottomValue);
                jacobi_row.push_back(i + shift);
                jacobi_column.push_back(cell.GetBottomIndex() + shift);

                sum = sum + bottomValue;
            }

            //left

            if (cell.Left == true){

                Cell &leftCell = cells[cell.GetLeftIndex()];

                double leftValue = 
                T_ij(cell.GetKx(), leftCell.GetKx(), HX) * 
                IfFuncRp(1 - cell.GetS(), 1 - leftCell.GetS(), cell.GetP(), leftCell.GetP()) + 
                K_si(PB, cell.GetP(), 1 - cell.GetS(), true) * 
                WI(cell.GetKx(), cell.GetKy());

                jacobi_value.push_back(-leftValue);
                jacobi_row.push_back(i + shift);
                jacobi_column.push_back(cell.GetLeftIndex() + shift);

                sum = sum + leftValue;
            }

            //right

             if (cell.Right == true){

                Cell &rightCell = cells[cell.GetRightIndex()];

                double rightValue = 
                T_ij(cell.GetKx(), rightCell.GetKx(), HX) * 
                IfFuncRp(1 - cell.GetS(), 1 - rightCell.GetS(), cell.GetP(), rightCell.GetP()) + 
                K_si(PB, cell.GetP(), 1 - cell.GetS(), true) * 
                WI(cell.GetKx(), cell.GetKy());

                jacobi_value.push_back(-rightValue);
                jacobi_row.push_back(i + shift);
                jacobi_column.push_back(cell.GetRightIndex() + shift);

                sum = sum + rightValue;
            }

            //cell 

            jacobi_value.push_back(sum);
            jacobi_row.push_back(i + shift);
            jacobi_column.push_back(i + shift);

        }
    }

    return 0;

}

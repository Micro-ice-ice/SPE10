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

    vector<double> jacobi_value; //COO jacobi matrix format
    vector<int> jacobi_row;
    vector<int> jacobi_column;

    //open data file

    string filename = "../data.txt";   // Name of the file
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
        //  -------------------
        //  |        |        |
        //  | dR0/dP | dR0/dS |
        //  |        |        |
        //  -------------------
        //  |        |        |
        //  | dRw/dS | dRw/dP |
        //  |        |        |
        //  -------------------
        //first block oil dR/dP
        {
            double sum = 0;

            double sj_top, sj_down, sj_left, sj_right;
            double pj_top, pj_down, pj_left, pj_right;

            //top 

            if (cell.Top){

                Cell &topCell = cells[cell.GetTopIndex()];

                double Tj = T_ij(cell.GetKy(), topCell.GetKy(), HY);

                double topValue = 
                Tj * IfFuncRp(cell.GetS(), topCell.GetS(), cell.GetP(), topCell.GetP()) +
                K_p(PB, cell.GetP(), cell.GetS()) *
                WI(cell.GetKx(), cell.GetKy());

                jacobi_value.push_back(-topValue);
                jacobi_row.push_back(i);
                jacobi_column.push_back(cell.GetTopIndex());

                sum = sum + Tj;
                pj_down = cell.GetP();
                sj_down = cell.GetS();
            }

            //bottom

            if (cell.Bottom){

                Cell &bottomCell = cells[cell.GetBottomIndex()];

                double Tj = T_ij(cell.GetKy(), bottomCell.GetKy(), HY);

                double bottomValue = 
                Tj * IfFuncRp(cell.GetS(), bottomCell.GetS(), cell.GetP(), bottomCell.GetP()) +
                K_p(PB, cell.GetP(), cell.GetS()) *
                WI(cell.GetKx(), cell.GetKy());

                jacobi_value.push_back(-bottomValue);
                jacobi_row.push_back(i);
                jacobi_column.push_back(cell.GetBottomIndex());

                sum = sum + Tj;
                pj_top = cell.GetP();
                sj_top = cell.GetS();
            }

            //left

            if (cell.Left){

                Cell &leftCell = cells[cell.GetLeftIndex()];

                double Tj = T_ij(cell.GetKx(), leftCell.GetKx(), HX);

                double leftValue = 
                Tj * IfFuncRp(cell.GetS(), leftCell.GetS(), cell.GetP(), leftCell.GetP()) +
                K_p(PB, cell.GetP(), cell.GetS()) *
                WI(cell.GetKx(), cell.GetKy());

                jacobi_value.push_back(-leftValue);
                jacobi_row.push_back(i);
                jacobi_column.push_back(cell.GetLeftIndex());

                sum = sum + Tj;
                pj_right = cell.GetP();
                sj_right = cell.GetS();
            }

            //right

             if (cell.Right){

                Cell &rightCell = cells[cell.GetRightIndex()];

                double Tj = T_ij(cell.GetKx(), rightCell.GetKx(), HX);

                double rightValue =
                Tj * IfFuncRp(cell.GetS(), rightCell.GetS(), cell.GetP(), rightCell.GetP()) +
                K_p(PB, cell.GetP(), cell.GetS()) *
                WI(cell.GetKx(), cell.GetKy());

                jacobi_value.push_back(-rightValue);
                jacobi_row.push_back(i);
                jacobi_column.push_back(cell.GetRightIndex());

                sum = sum + Tj;
                pj_left = cell.GetP();
                sj_left = cell.GetS();
            }

            //cell 

            jacobi_value.push_back(sum * IfFuncRpi(cell.GetS(), sj_top, sj_down, sj_left, sj_right,
                                                   cell.GetP(), pj_top, pj_down, pj_left, pj_right) + K_p(PB, cell.GetP(), cell.GetS())
                                                   * WI(cell.GetKx(), cell.GetKy()));
            jacobi_row.push_back(i);
            jacobi_column.push_back(i);

        }

        //second block oil dR/dS
        {
            int shift = NX * NY;
            double sum = 0;

            double sj_top, sj_down, sj_left, sj_right;
            double pj_top, pj_down, pj_left, pj_right;

            //top

            if (cell.Top){

                Cell &topCell = cells[cell.GetTopIndex()];

                double Tj = T_ij(cell.GetKy(), topCell.GetKy(), HY);

                double topValue =
                        Tj * (cell.GetP() - topCell.GetP()) * IfFuncRs(cell.GetP(), topCell.GetP(), false);

                jacobi_value.push_back(-topValue);
                jacobi_row.push_back(i);
                jacobi_column.push_back(cell.GetTopIndex());

                sum = sum + Tj;
                pj_down = cell.GetP();
                sj_down = cell.GetS();
            }

            //bottom

            if (cell.Bottom){

                Cell &bottomCell = cells[cell.GetBottomIndex()];

                double Tj = T_ij(cell.GetKy(), bottomCell.GetKy(), HY);

                double bottomValue =
                        Tj * (cell.GetP() - bottomCell.GetP()) * IfFuncRs(cell.GetP(), bottomCell.GetP(), false);

                jacobi_value.push_back(-bottomValue);
                jacobi_row.push_back(i);
                jacobi_column.push_back(cell.GetBottomIndex());

                sum = sum + Tj;
                pj_top = cell.GetP();
                sj_top = cell.GetS();
            }

            //left

            if (cell.Left){

                Cell &leftCell = cells[cell.GetLeftIndex()];

                double Tj = T_ij(cell.GetKx(), leftCell.GetKx(), HX);

                double leftValue =
                        Tj * (cell.GetP() - leftCell.GetP()) * IfFuncRs(cell.GetP(), leftCell.GetP(), false);

                jacobi_value.push_back(-leftValue);
                jacobi_row.push_back(i);
                jacobi_column.push_back(cell.GetLeftIndex());

                sum = sum + Tj;
                pj_right = cell.GetP();
                sj_right = cell.GetS();
            }

            //right

            if (cell.Right){

                Cell &rightCell = cells[cell.GetRightIndex()];

                double Tj = T_ij(cell.GetKx(), rightCell.GetKx(), HX);

                double rightValue =
                        Tj * (cell.GetP() - rightCell.GetP()) * IfFuncRs(cell.GetP(), rightCell.GetP(), false);

                jacobi_value.push_back(-rightValue);
                jacobi_row.push_back(i);
                jacobi_column.push_back(cell.GetRightIndex());

                sum = sum + Tj;
                pj_left = cell.GetP();
                sj_left = cell.GetS();
            }

            //cell

            double chosen_pj;
            double chosen_sj = IfFuncRpi(cell.GetS(), sj_top, sj_down, sj_left, sj_right,
                                  cell.GetP(), pj_top, pj_down, pj_left, pj_right);
            if (chosen_sj == sj_top){
                chosen_pj = pj_top;
            }
            else if (chosen_sj == sj_down){
                chosen_pj = pj_down;
            }
            else if (chosen_sj == sj_left){
                chosen_pj = pj_left;
            }
            else {
                chosen_pj = pj_right;
            }

            jacobi_value.push_back(cell.GetPhi() / HT + sum * (cell.GetP() - chosen_pj) * IfFuncRs(cell.GetP(), chosen_pj, true) +
                                                   K_s(PB, cell.GetP(), cell.GetS(), true)
                                                   * WI(cell.GetKx(), cell.GetKy() * (cell.GetP() - PB)));
            jacobi_row.push_back(i);
            jacobi_column.push_back(i + shift);

        }

        //third block water dR/dS
        {
            int shift = NX * NY;
            double sum = 0;

            double sj_top, sj_down, sj_left, sj_right;
            double pj_top, pj_down, pj_left, pj_right;

            //top

            if (cell.Top){

                Cell &topCell = cells[cell.GetTopIndex()];

                double Tj = T_ij(cell.GetKy(), topCell.GetKy(), HY);

                double topValue =
                        -Tj * (cell.GetP() - topCell.GetP()) * IfFuncRs(cell.GetP(), topCell.GetP(), false);

                jacobi_value.push_back(-topValue);
                jacobi_row.push_back(i);
                jacobi_column.push_back(cell.GetTopIndex());

                sum = sum + Tj;
                pj_down = cell.GetP();
                sj_down = cell.GetS();
            }

            //bottom

            if (cell.Bottom){

                Cell &bottomCell = cells[cell.GetBottomIndex()];

                double Tj = T_ij(cell.GetKy(), bottomCell.GetKy(), HY);

                double bottomValue =
                        -Tj * (cell.GetP() - bottomCell.GetP()) * IfFuncRs(cell.GetP(), bottomCell.GetP(), false);

                jacobi_value.push_back(-bottomValue);
                jacobi_row.push_back(i);
                jacobi_column.push_back(cell.GetBottomIndex());

                sum = sum + Tj;
                pj_top = cell.GetP();
                sj_top = cell.GetS();
            }

            //left

            if (cell.Left){

                Cell &leftCell = cells[cell.GetLeftIndex()];

                double Tj = T_ij(cell.GetKx(), leftCell.GetKx(), HX);

                double leftValue =
                        -Tj * (cell.GetP() - leftCell.GetP()) * IfFuncRs(cell.GetP(), leftCell.GetP(), false);

                jacobi_value.push_back(-leftValue);
                jacobi_row.push_back(i);
                jacobi_column.push_back(cell.GetLeftIndex());

                sum = sum + Tj;
                pj_right = cell.GetP();
                sj_right = cell.GetS();
            }

            //right

            if (cell.Right){

                Cell &rightCell = cells[cell.GetRightIndex()];

                double Tj = T_ij(cell.GetKx(), rightCell.GetKx(), HX);

                double rightValue =
                        -Tj * (cell.GetP() - rightCell.GetP()) * IfFuncRs(cell.GetP(), rightCell.GetP(), false);

                jacobi_value.push_back(-rightValue);
                jacobi_row.push_back(i);
                jacobi_column.push_back(cell.GetRightIndex());

                sum = sum + Tj;
                pj_left = cell.GetP();
                sj_left = cell.GetS();
            }

            //cell

            double chosen_pj;
            double chosen_sj = IfFuncRpi(cell.GetS(), sj_top, sj_down, sj_left, sj_right,
                                         cell.GetP(), pj_top, pj_down, pj_left, pj_right);
            if (chosen_sj == sj_top){
                chosen_pj = pj_top;
            }
            else if (chosen_sj == sj_down){
                chosen_pj = pj_down;
            }
            else if (chosen_sj == sj_left){
                chosen_pj = pj_left;
            }
            else {
                chosen_pj = pj_right;
            }

            jacobi_value.push_back(-cell.GetPhi() / HT - sum * (cell.GetP() - chosen_pj) * IfFuncRs(cell.GetP(), chosen_pj, true) +
                                   K_s(PB, cell.GetP(), cell.GetS(), true)
                                   * WI(cell.GetKx(), cell.GetKy() * (PB - cell.GetP())));
            jacobi_row.push_back(i + shift);
            jacobi_column.push_back(i);

        }

        //forth block water dR/dP
        {
            int shift = NX * NY;
            double sum = 0;

            double sj_top, sj_down, sj_left, sj_right;
            double pj_top, pj_down, pj_left, pj_right;

            //top

            if (cell.Top){

                Cell &topCell = cells[cell.GetTopIndex()];

                double Tj = T_ij(cell.GetKy(), topCell.GetKy(), HY);

                double topValue =
                        Tj * IfFuncRp(1 - cell.GetS(), 1 - topCell.GetS(), cell.GetP(), topCell.GetP()) +
                        K_p(PB, cell.GetP(), cell.GetS()) *
                        WI(cell.GetKx(), cell.GetKy());

                jacobi_value.push_back(-topValue);
                jacobi_row.push_back(i + shift);
                jacobi_column.push_back(cell.GetTopIndex() + shift);

                sum = sum + Tj;
                pj_down = cell.GetP();
                sj_down = cell.GetS();
            }

            //bottom

            if (cell.Bottom){

                Cell &bottomCell = cells[cell.GetBottomIndex()];

                double Tj = T_ij(cell.GetKy(), bottomCell.GetKy(), HY);

                double bottomValue =
                        Tj * IfFuncRp(1 - cell.GetS(), 1 - bottomCell.GetS(), cell.GetP(), bottomCell.GetP()) +
                        K_p(PB, cell.GetP(), cell.GetS()) *
                        WI(cell.GetKx(), cell.GetKy());

                jacobi_value.push_back(-bottomValue);
                jacobi_row.push_back(i + shift);
                jacobi_column.push_back(cell.GetBottomIndex() + shift);

                sum = sum + Tj;
                pj_top = cell.GetP();
                sj_top = cell.GetS();
            }

            //left

            if (cell.Left){

                Cell &leftCell = cells[cell.GetLeftIndex()];

                double Tj = T_ij(cell.GetKx(), leftCell.GetKx(), HX);

                double leftValue =
                        Tj * IfFuncRp(1 - cell.GetS(), 1 - leftCell.GetS(), cell.GetP(), leftCell.GetP()) +
                        K_p(PB, cell.GetP(), cell.GetS()) *
                        WI(cell.GetKx(), cell.GetKy());

                jacobi_value.push_back(-leftValue);
                jacobi_row.push_back(i + shift);
                jacobi_column.push_back(cell.GetLeftIndex() + shift);

                sum = sum + Tj;
                pj_right = cell.GetP();
                sj_right = cell.GetS();
            }

            //right

            if (cell.Right){

                Cell &rightCell = cells[cell.GetRightIndex()];

                double Tj = T_ij(cell.GetKx(), rightCell.GetKx(), HX);

                double rightValue =
                        Tj * IfFuncRp(1 - cell.GetS(), 1 - rightCell.GetS(), cell.GetP(), rightCell.GetP()) +
                        K_p(PB, cell.GetP(), cell.GetS()) *
                        WI(cell.GetKx(), cell.GetKy());

                jacobi_value.push_back(-rightValue);
                jacobi_row.push_back(i + shift);
                jacobi_column.push_back(cell.GetRightIndex() + shift);

                sum = sum + Tj;
                pj_left = cell.GetP();
                sj_left = cell.GetS();
            }

            //cell

            jacobi_value.push_back(sum * IfFuncRpi(1 - cell.GetS(), 1 - sj_top, 1 - sj_down, 1 - sj_left, 1 - sj_right,
                                                   cell.GetP(), pj_top, pj_down, pj_left, pj_right) +
                                                   K_s(PB, cell.GetP(), cell.GetS(), true) * WI(cell.GetKx(), cell.GetKy()));
            jacobi_row.push_back(i + shift);
            jacobi_column.push_back(i + shift);

        }

    }

    cout << jacobi_value.size();

    return 0;

}

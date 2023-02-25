#pragma once

class Cell
{

private:

    static int Indexator;

    static int Nx;

    static int Ny;

    static int CellCount;

    int Index;

    int TopIndex;

    int BottomIndex;

    int LeftIndex;

    int RightIndex;
    
    double Kx;

    double Ky;

    double Phi;

public:

    bool Top;

    bool Bottom;

    bool Left;

    bool Right;

    static void SetNx(int nx){

        Nx = nx;
    }

    static void SetNy(int ny){

        Ny = ny;
    }

    static void SetCellCount(int count){

        CellCount = count;
    }

    int GetTopIndex(){

        return TopIndex;
    }

    int GetBottomIndex(){

        return BottomIndex;
    }

    int GetLeftIndex(){

        return LeftIndex;
    }

    int GetRightIndex(){

        return RightIndex;
    }

    // void SetTopIndex(int index){

    //     TopIndex = index;
    // }

    // void SetBottomIndex(int index){

    //     BottomIndex = index;
    // }

    // void SetLeftIndex(int index){

    //     LeftIndex = index;
    // }

    // void SetRightIndex(int index){

    //     RightIndex = index;
    // }

    double GetKx(){

        return Kx;
    }

    double GetKy(){

        return Ky;
    }

    double GetPhi(){

        return Phi;
    }

    Cell(double kx, double ky, double phi);

    ~Cell();
};

Cell::Cell(double kx, double ky, double phi)
{
    Index = Indexator;
    Indexator++;
    Kx = kx;
    Ky = ky;
    Phi = phi;

    //top

    if (Index - Nx >= 0){

        TopIndex = Index - Nx;
        Top = true;
    }
    else{

        Top = false;
    }

    //bottom

    if (Index + Nx <= CellCount){

        BottomIndex = Index + Nx;
        Bottom = true;
    }
    else{

        Bottom = false;
    }

    //left

    if (Index % Nx != 0){

        LeftIndex = Index - 1;
        Left = true;
    }
    else{

        Left = false;
    }

    //right

    if ((Index + 1) % Nx != 0){

        RightIndex = Index + 1;
        Right = true;
    }
    else{

        Right = false;
    }

}

Cell::~Cell()
{
}

//init static vars

int Cell::Indexator = 0;
int Cell::Nx = 0;
int Cell::Ny = 0;
int Cell::CellCount = 0;
#ifndef SSATOGGLE_H_

#define SSATOGGLE_H_
#include "Cells.h"

struct ToggleCell{
    int lineage;
    int parent;
    int rcdSize;
    double* time;
    int* green;
    int* red;
    // Initializes the toggle state of a cell, and cell states are blank.
    ToggleCell(){
        this->lineage = -1;
        this->parent = -1;
        this->rcdSize = 0;
        this->time = nullptr;
        this->green = nullptr;
        this->red = nullptr;
    }
    // Initializes the Toggle state of a cell, this function is identical to initCell(...).
    ToggleCell(const double& startTime, const double &endTime, const double &outputTime,
               const int& initgreen, const int& initred, const int& parent, const int& lineage){
        int runsize = (int) floor((endTime-startTime ) / outputTime) +1;
        this->green = new int[runsize];
        this->red = new int[runsize];
        this->time = new double [runsize];
        *(this->green) = initgreen;
        *(this->red) = initred;
        *(this->time) = startTime;
        this->parent = parent;
        this->lineage= lineage;
        this->rcdSize = 1;
    }
    ~ToggleCell(){
        if(this->time){
            delete[] this->time;
            delete[] this->green;
            delete[] this->red;
            this->time = nullptr;
            this->green = nullptr;
            this->red = nullptr;
        }
    }
    // deep copy
    ToggleCell& operator=(ToggleCell& cellTemp){
        lineage = cellTemp.lineage;
        parent = cellTemp.parent;
        rcdSize = cellTemp.rcdSize;
        if(cellTemp.rcdSize > 0){
            time = new double[rcdSize];
            green = new int[rcdSize];
            red = new int[rcdSize];
            for(int i=0; i<cellTemp.rcdSize; ++i){
                green[i] = cellTemp.green[i];
                red[i] = cellTemp.red[i];
                time[i] = cellTemp.time[i];
            }
        }
        return *this;
    };

};

struct cellBatch{
    ToggleCell* cells;
    int size;
    cellBatch(){
        this->cells = nullptr;
        this->size = 0;
    }
    ~cellBatch(){
        if(this->cells){
            delete[] this->cells;
            this->cells = nullptr;
        }
    }

};

//int runSim(const double& gr, const double& endTime, double* green, double* red);
//int runSim(const double& gr, const int& green, const int& red,
//           const double& endTime, const double& outputTime, const double& t0,
//           double* saveT, int* saveX1, int* saveX2, int* saveSize);
int runSim(const double& timeNow, Cell* pCell);
int rumMultiSim(const int& threadNum, const double& gr, int* green, int* red,
                const double& endTime, const double& outputTime, int simsize, double* saveBuff, int* saveLength);
void runBatchSim(const int& threadNum, const double& gr, const int& green, const int& red,
                 const double& endTime, const double& outputTime, const int &maxCell,
                 ToggleCell** cellsarray, int* cellsSize);
void freeCellArray(ToggleCell* cell, int& size);
void freeCellMem(ToggleCell* cell);


double hillFunc(const double &leakage, const double &k, const double &n, const double &p);
double alphaG(const double &gr);
double alphaR(const double &gr);

#endif  /* End SSATOGGLE_H_ */
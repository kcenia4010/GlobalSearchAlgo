#include "wrappers_functions.h"

void CreateGrishaginFunction(){
    myClassGrishaginFunction = new vagrish::GrishaginFunction();
}

void MySetFunctionNumber(int i){
    myClassGrishaginFunction->SetFunctionNumber(i);
}

double GetBoundsA0(){
    double lb[2] = {0 ,0};
    double ub[2] = {0 ,0};
    myClassGrishaginFunction->GetBounds(lb, ub);
    return lb[0];
}

double GetBoundsA1(){
    double lb[2] = {0 ,0};
    double ub[2] = {0 ,0};
    myClassGrishaginFunction->GetBounds(lb, ub);
    return lb[1];
}

double GetBoundsB0() {
    double lb[2] = {0 ,0};
    double ub[2] = {0 ,0};
    myClassGrishaginFunction->GetBounds(lb, ub);
    return ub[0];
}

double GetBoundsB1() {
    double lb[2] = {0 ,0};
    double ub[2] = {0 ,0};
    myClassGrishaginFunction->GetBounds(lb, ub);
    return ub[1];
}

double MyCalculate(double y0, double y1){
    double y[2] = {y0, y1};
    return myClassGrishaginFunction->Calculate(y);
}

void FreeGrishaginFunction(){
    myClassGrishaginFunction->~GrishaginFunction();
}
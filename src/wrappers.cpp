#include "wrappers.h"

void CreateEvovlent() {
    myClassEvovlent = new TEvolvent(2, 10);
}

void SetBounds(double A0, double A1, double B0, double B1){
    double _A[2];
    double _B[2];
    _A[0] = A0; _A[1] = A1; 
    _B[0] = B0; _B[1] = B1; 
    myClassEvovlent->SetBounds(_A, _B);
}

double GetImage0(double x){
    myClassEvovlent->GetImage(x, result_image);
    return result_image[0];
}

double GetImage1(double x){
    return result_image[1];
}

void FreeEvolvent(){
    myClassEvovlent->~TEvolvent();
}
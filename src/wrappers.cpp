#include "wrappers.h"

void CreateEvovlent(int idx) {
    threads[idx].myClassEvovlent = new TEvolvent(2, 10);
}

void CopyEvovlent(int idx) {
    threads[idx].myClassEvovlent = new TEvolvent(*threads[0].myClassEvovlent);
}

void SetBounds(int idx, double A0, double A1, double B0, double B1){
    double _A[2];
    double _B[2];
    _A[0] = A0; _A[1] = A1; 
    _B[0] = B0; _B[1] = B1; 
    threads[idx].myClassEvovlent->SetBounds(_A, _B);
}

double GetImage0(int idx, double x){
    threads[idx].myClassEvovlent->GetImage(x, threads[idx].result_image);
    return threads[idx].result_image[0];
}

double GetImage1(int idx, double x){
    return threads[idx].result_image[1];
}

void FreeEvolvent(int idx){
    threads[idx].myClassEvovlent->~TEvolvent();
}
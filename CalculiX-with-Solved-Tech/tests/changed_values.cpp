// A correctness check of two coefficient updates, never a performance case.
#include <cmath>
#include <cstdio>
extern "C" void spooles(double*,double*,double*,double*,double*,double*,int*,int*,int*,int*,int*,int*,int*);
int main(){
    int n=3,nz=2,zero=0,column_counts[]={1,1,0},row_ids[]={2,3};double shift=0;
    double diagonal[]={4,5,6},off_diagonal[]={1,2},rhs[]={6,17,22};
    spooles(diagonal,off_diagonal,nullptr,nullptr,&shift,rhs,column_counts,row_ids,&n,&nz,&zero,&zero,&nz);
    for(int i=0;i<n;++i)if(std::abs(rhs[i]-(i+1))>1e-10)return 91;
    for(double &x:diagonal)x*=2;
    double updated_rhs[]={10,27,40};
    spooles(diagonal,off_diagonal,nullptr,nullptr,&shift,updated_rhs,column_counts,row_ids,&n,&nz,&zero,&zero,&nz);
    for(int i=0;i<n;++i)if(std::abs(updated_rhs[i]-(i+1))>1e-10)return 92;
    puts("CHANGED_VALUES_PASSED: same graph, changed coefficients/RHS, known solution (1,2,3)");
}

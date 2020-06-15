
#include <Eigen/Core>
#include <iostream>

#include "ddgsolver/util.h"

int main() {
  Eigen::Array<double, 4, 3> A, B;

  for(int i = 0; i < 4; ++i){
    A(i,0) = 3*i;
    A(i,1) = 3*i+1;
    A(i,2) = 3*i+2;
  }

  for(int i = 0; i < 4; ++i){
    B(i,0) = 10*i;
    B(i,1) = 10*i+1;
    B(i,2) = 10*i+2;
  }
  std::cout << "A:" << std::endl << A << std::endl;
  std::cout << "B:" << std::endl << B << std::endl; 

  std::cout << "Sizeof(A): " << A.rows() << "x" << A.cols() << std::endl;
  std::cout << "Sizeof(B): " << B.rows() << "x" << B.cols() << std::endl;
  std::cout << A + B << std::endl;

  auto res = ddgsolver::dot(A,B);
  std::cout << "Sizeof(res): " << res.rows() << "x" << res.cols() << std::endl;
  std::cout << res << std::endl; 

  auto res2 = (A.cwiseProduct(B).rowwise().sum());
  std::cout << "Sizeof(res2): " << res2.rows() << "x" << res2.cols() << std::endl;
  std::cout << res2 << std::endl;
  return 0;
}

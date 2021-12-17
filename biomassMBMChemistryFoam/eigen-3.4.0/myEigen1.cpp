#include <iostream>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/IterativeLinearSolvers>



int main()
{
  int n = 10;
  Eigen::SparseMatrix<double> S = Eigen::MatrixXd::Random(n,n).sparseView(0.5,1);
  S = S.transpose()*S;
  Eigen::VectorXd b(n), x;
  b.setRandom();
  Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Upper> cg;

  x = cg.compute(S).solve(b);
  std::cout << "CG:       #iterations: " << cg.iterations() << ", estimated error: " << cg.error() << std::endl;

  // Eigen::SparseMatrix<double> A;
  // Eigen::MatrixXd m(2,2);
  // Eigen::VectorXd b, x;
  // Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Upper> solver;
  // x = solver.compute(A).solve(b);

  // m(0,0) = 3;
  // m(1,0) = 2.5;
  // m(0,1) = -1;
  // m(1,1) = m(1,0) + m(0,1);
  std::cout << x << std::endl;
}

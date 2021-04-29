#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <dune/common/parallel/mpihelper.hh>

#include <dune/istl/bcrsmatrix.hh> 
#include <dune/istl/bvector.hh>    
#include <dune/istl/operators.hh>
#include <dune/istl/preconditioners.hh>
#include <dune/istl/solvers.hh>

#include <dune/grid/io/file/vtk/vtksequencewriter.hh>
#include <dune/grid/io/file/vtk.hh>
#include <dune/grid/onedgrid.hh>

#include <cassert>
#include <cmath>
#include <iostream>
#include <string>

/// Duljina intervala
const double L = 2.0;
/// Difuzijski koeficijent
const double D = 0.001;
/// Vrijeme simulacije
const double T = 10.0;
/// Broj elemenata u mreži
static unsigned int n_elements = 50;
/// Desna strana.
double f(double x, double t) {
    return  (5*L*x)/(2*M_PI)*cos(2*M_PI*x*t/L) + 5*D*sin(2*M_PI*x*t/L);
}
/// Egzaktno rješenje
double sol(double x, double t) {
    return (5*L*L)/(4*M_PI*M_PI)*sin(2*M_PI*x*t/L) + 2*x/L;
}

int main(int argc, char **argv) {
    Dune::MPIHelper::instance(argc, argv);

    using Scalar = Dune::FieldVector<double, 1>;
    using Vector = Dune::BlockVector<Scalar>;
    using Matrix = Dune::BCRSMatrix<Dune::FieldMatrix<double, 1, 1>>;

    using GridType = Dune::OneDGrid;
    using GridView = GridType::LeafGridView;

    if (argc > 1){
        n_elements = std::stoi(argv[1]);
    }
    assert(n_elements > 0);
    //AX=F
    Scalar D;
    Vector F, X, X_exact;
    Matrix A;

    GridType grid(n_elements, 0, L);
    GridView gridview = grid.leafGridView();

    auto n_vertices = gv.size(1);

    assert(n_elements + 1 = n_vertices);

    F.resize(n_vertices);
    X.resize(n_vertices);
    X_exact.resize(n_vertices);
    X = 0.0;

    double h = L/n_elements;
    double h2 = h * h;
  
    //rubovi??? svakim korakom je drugi delta t SVE NEKAKO STAVITI U WHILE PETLJU HMMMMMMMM
    //initialize matrix A
    A.setSize(n_vertices,n_vertices);
    A.setBuildMode(Matrix::random);
    A.setrowsize(0,1); //ili?
    for(unsigned int i = 1; i < n_vertices-1; ++i){
        A.setrowsize(i,3);
    }
    A.setrowsize(n_vertices-1,1);
    A.endrowsizes();

    A.addindex(0,0);
    for(unsigned int i = 1; i < n_vertices-1; ++i){
        A.addindex(i,i-1);
        A.addindex(i,i);
        A.addindex(i,i+1);
    }
    A.addindex(n_vertices-1,n_vertices-1);
    A.endindices();
    double lambda = D * dt / h2;
    A[0][0] = 1.0;
    for(unsigned int i = 1; i < n_vertices-1; ++i){
        A[i][i-1] = -lambda;
        A[i][i] = 1 + 2*lambda;
        A[i][i+1] = -lambda;
    }
    A[n_vertices-1][n_vertices-1] = 1.0;
    //initialize vector F
    F[0] = ;
    F[n_vertices - 1] = ;
    for(unsigned int i = 0; i < n_vertices; ++i){
        F[i] = dt * f(i*h);
    }
    for(unsigned int i = 0; i < n_vertices;++i){
        X_exact[i] = sol(i*h);
    }
    Dune::MatrixAdapter<Matrix,Vector,Vector> op(A);
    Dune::SeqILU<Matrix,Vector,Vector> prec(A, 0.92);

    Dune::BiCGSTABSolver<Vector> solver(op, prec, 1e-8, 1000, 10);
    Dune::InverseOperatorResult r;

    auto FF = F;
    solver.apply(X,FF,r);

    std::cout << "Did solver converge?  " ;
    std::cout << r.converged << std::endl;
    if(r.converged){
      std::cout << "Number of iterations needed: " << std::endl;
      std::cout << r.reduction << std::endl;
    }

    Dune::VTKSequenceWriter<GridView> writer(gridview,"solution","x","t");
    writer.addVertexData(xnew, "numeric");
    writer.addVertexData(X_exact, "exact");
    writer.write("solution");

    return 0;
}
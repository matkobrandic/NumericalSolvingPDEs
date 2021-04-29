#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <dune/common/parallel/mpihelper.hh>

#include <dune/istl/bvector.hh>      // BlockVector
#include <dune/istl/bcrsmatrix.hh>   // Blok matrica
#include <dune/istl/operators.hh>
#include <dune/istl/preconditioners.hh>
#include <dune/istl/solvers.hh>

#include <dune/grid/io/file/vtk.hh>
#include <dune/grid/onedgrid.hh>

#include <cmath>
#include <cassert>
#include <string>

/**
 *  Rješavamo zadaću
 *  \[ -u'' = f \quad \text{na } (0,L)\]
 *  \[  u(o) = g_0,\; u(L)=g_1 \]
 *  metodom konačnih diferencija s uniformnim prostornim korakom h.
 */

/// Duljina intervala
const double L   = 2.0;
///  Rubni uvjet
const double g_0 = 0.0;
const double g_1 = 2.0;
/// Broj elemenata u mreži
static unsigned int n_elem = 10;

/// Desna strana.
double f(double x){
    return 5*std::sin(2*M_PI*x/L);
}
/// Egzaktno rješenje
double sol(double x){
    return 5*L*L*std::sin(2*M_PI*x/L)/(4*M_PI*M_PI) + (g_1-g_0)*x/L + g_0;
}

int main(int argc, char *argv[]){

    Dune::MPIHelper::instance(argc, argv);

    if(argc > 1){
        n_elem = std::stoi(argv[1]);
    }

    using Vector = Dune::BlockVector< Dune::FieldVector<double, 1> >;
    using Matrix = Dune::BCRSMatrix< Dune::FieldMatrix<double, 1, 1> >;
    //AX=F
    Vector F, X, X_exact;
    Matrix A;

    using GridType = Dune::OneDGrid;
    using GridView = GridType::LeafGridView;

    GridType grid(n_elem, 0, L);
    GridView gv = grid.leafGridView();

    auto n_vertices = gv.size(1);

    assert(n_elem + 1 = n_vertices);

    F.resize(n_vertices);
    X.resize(n_vertices);
    X_exact.resize(n_vertices);
    X = 0.0;

    double h = L/n_elem;
    double h2 = h * h;
    //Vector F
    F[0] = g_0;
    F[n_vertices-1] = g_1;
    for(unsigned int i = 1; i < n_vertices-1; ++i){
        F[i] = h2 * f(i*h);
    }
    for(unsigned int i = 0; i < n_vertices;++i){
        X_exact[i] = sol(i*h);
    }
    //Matrix A
    A.setSize(n_vertices,n_vertices);
    A.setBuildMode(Matrix::random);
    A.setrowsize(0,1);
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

    A[0][0] = 1.0;
    for(unsigned int i = 1; i < n_vertices-1; ++i){
        A[i][i-1] = -1;
        A[i][i] = 2;
        A[i][i+1] = -1;
    }
    A[n_vertices-1][n_vertices-1] = 1.0;

    Dune::MatrixAdapter<Matrix,Vector,Vector> op(A);
    Dune::SeqILU<Matrix,Vector,Vector> prec(A, 0.91); //odi u dune-istl dokumentaciju i nađi SeqILU

    Dune::BiCGSTABSolver<Vector> solver(op, prec, 1e-8,1000,10);
    Dune::InverseOperatorResult r;

    auto FF = F;
    solver.apply(X,FF,r);



    if(r.converged){
        std::cout << r.converged << std::endl;
        std::cout << r.iterations << std::endl;
        std::cout << r.reduction << std::endl;
    }
    else{
        std::cout << "Solver did not converge!\n" << std::endl;
    }

    Dune::VTKWriter<GridView> writer(gv);
    writer.addVertexData(X,"solved solution");
    writer.addVertexData(X_exact,"exact solution");
    writer.write("solution");


    return 0;
}

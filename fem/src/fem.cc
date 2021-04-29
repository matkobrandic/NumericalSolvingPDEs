#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <array>
#include <cmath>

#include <dune/geometry/quadraturerules.hh>

#include <dune/grid/io/file/dgfparser/dgfug.hh>
#include <dune/grid/io/file/dgfparser/gridptr.hh>
#include <dune/grid/io/file/gmshreader.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/grid/uggrid.hh>

#include <dune/common/fvector.hh>
#include <dune/grid/common/mcmgmapper.hh>

#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/bvector.hh>
#include <dune/istl/ilu.hh>
#include <dune/istl/io.hh>
#include <dune/istl/operators.hh>
#include <dune/istl/preconditioners.hh>
#include <dune/istl/solvers.hh>

#include <dune/localfunctions/lagrange/pk.hh>

/*
  Kod za rješavanje rubne zadaće:
    - Laplace u + a0(x)u = f(x),   x\in Omega
                       u = g      na \Gamma_D
            \nabla u . n = h      na \Gamma_N

  Kod će biti testiran na egzaktnom rješenju.
  Zadaća se riješava P1 elementima na simpleksima.
*/


// Egzaktno rješenje kao funkcija od globalne koordinate
// Koristimo egzaktno rješenje radi testiranja koda.
template <int dim>
double exact(Dune::FieldVector<double, dim> x){
//  double result = 1.0 + 3 * (x * x);
    double result = 10*std::sin(2*M_PI*x[0]);
    for(int i = 1; i < dim; ++i){
        result *= std::sin(2*M_PI*x[i]);
    }
    return result;
}
// Slobodni koeficijent
template <int dim>
double a0(Dune::FieldVector<double, dim> x){
    double result = 1.0 + x*x;
    return result;
}
// Desna strana jednadžbe.
// Računamo ju koristeći egzaktno rješenje pomoću konačnih diferencija.
template <int dim>
double f(Dune::FieldVector<double, dim> x){
    std::array<Dune::FieldVector<double, dim>, dim> xp;
    std::array<Dune::FieldVector<double, dim>, dim> xm;
    std::array<double, dim> cp;
    std::array<double, dim> cm;
    const double h = 1E-6;  // Korak diferenciranja
    const double h2 = h*h;
    for(unsigned int i = 0; i < dim; ++i){
        xp[i] = x;
        xm[i] = x;
        xp[i][i] += h;
        xm[i][i] -= h;
        cp[i] = exact(xp[i]);
        cm[i] = exact(xm[i]);
    }
    double c0 = exact(x);

    // operator - Laplace aproksimiran konačnim diferencijama u točki x.
    double result = 0.0;
    for (int i = 0; i < dim; i++){
        result -= (cp[i] - 2 * c0 + cm[i]) / h2;
    }
    // dodaj slobodni član
    result += a0(x) * c0;

    return result;
}

// Identificiraj Neumannovu granicu
template <int dim>
bool isNeumann(Dune::FieldVector<double, dim> x){
    assert(dim == 3);
    const double EPS = 1E-10;
    //    Neumann za z=0 i z=1
    if(std::abs(x[2]) < EPS or std::abs(x[2] - 1.0) < EPS){
        return true;
    }
    return false;
}

// Identificiraj Dirichletovu granicu
template <int dim>
bool isDirichlet(Dune::FieldVector<double, dim> x){
    return !isNeumann(x);
}

// Funkcija vraća Neumannov rubni uvjet.
// x = toka na Neumannovoj granici
// unormal = jedinična vanjska normala na granicu
template <int dim>
double flux(Dune::FieldVector<double, dim> x,
            Dune::FieldVector<double, dim> unormal){
    const double h = 1E-6;  // Korak diferenciranja
    Dune::FieldVector<double, dim> x1 = unormal;
    x1 *= h;
    x1 += x;
    double u1 = exact(x1);
    double u = exact(x);
    return (u1 - u) / h;
}
// GV = GridView
// Vec = vektor (=BlockVector<...>)
// Izračunaj egzaktno rješenje u vrhovima mreže i vrati izračunate vrijednosti
// kroz coeff.
template <class GV, class Vec>
void vertexdata(const GV &gridView, Vec &coeff){
    const int dim = GV::dimensionworld;
    auto const & idset = gridView.indexSet();
    coeff.resize(idset.size(dim));
    for(auto const & vertex : vertices(gridView)){
        const auto global = vertex.geometry().corner(0);
        const auto idx = idset.index(vertex);
        coeff[idx] = exact(global);
    }
}
// Izračunaj  matrixProfile. Ovdje ćemo radi preglednosti koda prvo naći
// indekse svih susjeda elementa. Pri tome koristimo skup std::set.
// Prednost skupa je što ignorira duplikate.
//  matrixProfile[i] = skup svih indeksa vrhova koji su susjedni vrhu i.
//                     Biti susjedan znači pripadati istom elementu mreže.
// GV = GridView
template <class GV>
void computeMatrixProfile(GV const &gv, std::vector<std::set<int>> &matrixProfile){
    const int dim = GV::dimensionworld;
    // Broj vrhova mreže
    const auto N = gv.size(dim);
    matrixProfile.resize(N);
    // Mapper koji indeksira sve vrhove mreže
    auto const &idset = gv.indexSet();

    // Po svim elementima mreže
    for(auto const &element : elements(gv)){
        // Broj vrhova u elementu
        int verticesNumber = element.geometry().corners();
        // A_{i,j} je različito od nule samo ako vrhovi s indeksima i,j pripadaju
        // istom elementu. Stoga je dovoljno na svakom element skupiti sve parove
        // indeksa vrhova.
        for(int i = 0; i < verticesNumber; ++i){
            auto index_i = idset.subIndex(element, i, dim);
            for(int j = 0; j < verticesNumber; ++j){
                auto index_j = idset.subIndex(element, j, dim);
                matrixProfile[index_i].insert(index_j);
            }
        }
    }
    std::cout << " Profil matrice je izračunat.\n";
}
// Dimenzioniraj i inicijaliziraj nulama matricu krutosti i vektor
// desne strane. Matrica krutosti se dimenzionira prema podacima u
// varijabli matrixProfile. Metoda  izracunajmatrixProfile() mora biti
// pozvana prije ove metode pa je pozivamo ovdje.
template <typename GV, typename Mat, typename Vec>
void initialize(GV const &gv, Mat &A, Vec &b){
    const int dim = GV::dimensionworld;
    std::vector<std::set<int>> matrixProfile;
    computeMatrixProfile(gv, matrixProfile);
    // Broj vrhova mreže
    const auto N = gv.size(dim);
    // dimenzije matrice A i vektora desne strane b
    A.setSize(N, N);
    A.setBuildMode(Mat::random);
    b.resize(N, false);
    // Postavimo broj ne-nul elemenata u svakom retku
    for(int i = 0; i < N; ++i){
        A.setrowsize(i, matrixProfile[i].size());
    }
    A.endrowsizes();  // dimenzioniranje redaka kompletirano
    // Definiraj indekse stupaca u svakom retki
    for(int i = 0; i < N; ++i){
        for(auto j : matrixProfile[i]){
            A.addindex(i, j);
        }
    }
    A.endindices();  // definiranje stupaca završeno
    // Profil matrice je time fiksiran
    // inicijalizacija nulama
    A = 0.0;
    b = 0.0;
    std::cout << " Matrica i vektor desne strane su inicijalizirani.\n";
}
// Asembliranje matrice krutosti i vektora desne strane.
// Profil matrice mora biti definiran i matrica i vektor desne strane moraju
// inicijalizirani nulama. To radi metoda init() koja se poziva na početku.
// GV = LeafGridView
// FEM = Lokalni konačni element
// Mat = tip matrice
// Vec = tip vektora
template <typename GV, typename FEM, typename Mat, typename Vec>
void assemble(GV &gv, FEM const &fem, Mat &A, Vec &b){
    // A = matrica krutosti, b = vektor desne strane
    const int dim = GV::dimensionworld;
    using FEMGradient = typename FEM::Traits::LocalBasisType::Traits::JacobianType;
    using FEMRange    = typename FEM::Traits::LocalBasisType::Traits::RangeType;
    initialize(gv, A, b);
    auto const & mapper = gv.indexSet();
    //asembliranje volumnog doprinosa
    for(auto const & element : elements(gv)){
        auto gt = element.type();
        auto vertexSize = element.geometry().corners();
        auto basisSize = fem.localBasis().size();
        // provjera da baznih funkcija ima koliko i vrhova
        assert(vertexSize == basisSize);
        // kvadraturna formula
        const auto &rule = Dune::QuadratureRules<double, dim>::rule(gt, 2);
        for(auto const & qpoint : rule){
            auto const weight = qpoint.weight();
            auto const & point = qpoint.position();
            auto dx = element.geometry().integrationElement(point);
            //Determinanta Jakobijana. integrationElement uzima lokalnu točku koju želimo integrirati
            auto jacobian = element. geometry().jacobianInverseTransposed(point);
            auto acoeff = a0(element.geometry().global(point));
            auto fcoeff = f(element.geometry().global(point));
            //sve bazne funkcije i gradijenti računaju se odjednom -> izvan petlji
            std::vector<FEMRange> phiHat(basisSize);
            fem.localBasis().evaluateFunction(point, phiHat);

            std::vector<FEMGradient> gradPhiHat(basisSize);
            fem.localBasis().evaluateJacobian(point, gradPhiHat);
            //  Po svim baznim funkcijama
            for(unsigned int i = 0; i < basisSize; i++){
                auto index_i = mapper.subIndex(element, i, dim);

                b[index_i] += fcoeff * phiHat[i] * dx * weight;  // volumni dio desne strane
                Dune::FieldVector<double, dim> gradient1;
                jacobian.mv(gradPhiHat[i][0], gradient1); //gradient = jacobian*gradPhiHat[i][0]
                // Po svim baznim funkcijama
                for(int j = 0; j < basisSize; j++){
                    auto index_j = mapper.subIndex(element, j, dim);
                    Dune::FieldVector<double, dim> gradient2;
                    jacobian.mv(gradPhiHat[j][0], gradient2);
                    A[index_i][index_j] += ((gradient1*gradient2) +
                                            acoeff*phiHat[i]*phiHat[j]) * dx * weight;
                }
            }
        }
    } // Kraj petlje po elementima
    // Dirichletov rubni uvjet
    // Svaku Dirichletovu točku zamijenjujemo trivijalnim linijama.
    for(auto const &element : elements(gv)){
        // instanciramo referentni element
        Dune::GeometryType gt = element.type();
        auto const & ref = Dune::referenceElement<double, dim>(gt);
        auto basisSize = ref.size(dim);  // broj vrhova = broj baznih funkcija

        for(auto const & side : intersections(gv, element)){
            // jesmo li na granici
            if(side.boundary()){
                // Broj vrhova na stranici
                //indexInInside gleda je li stranica na rubu ili ne
                //isto tako i sa indexInOutside
                //ref.size(f,codim,dim) - f,codim određuju subentitet
                //dakle nama f,codim daju u ovom slučaju stranicu
                //koliko ta stranica ima elemenata veličine dim??
                auto verticesNumber = ref.size(side.indexInInside(), 1, dim);
                assert(verticesNumber == side.geometry().corners());
                auto center = side.geometry().center();  // centar stranice
                if(isNeumann(center)){
                    auto gt = side.geometry().type();
                    const auto & rule = Dune::QuadratureRules<double, dim-1>::rule(gt, 2);
                    const auto unitNormal = side.centerUnitOuterNormal();
                    for(auto const & qpoint : rule){
                        //isti je kod, ali se nalazimo na
                        //referentnom elementu stranice, a ne elementa (dim-1)
                        auto const weight = qpoint.weight();
                        auto const & point = qpoint.position();
                        auto const dy = side.geometry().integrationElement(point);
                        const auto globalPoint = side.geometry().global(point);
                        const auto localPoint = element.geometry().local(globalPoint);
                        const auto hcoeff = flux(globalPoint, unitNormal);
                        std::vector<FEMRange> phiHat(basisSize);
                        fem.localBasis().evaluateFunction(localPoint, phiHat);
                        for(int i = 0; i < verticesNumber; i++){
                            //indexi sa stranicom na referentnom elementu - predavanje
                            auto local_index_i = ref.subEntity(side.indexInInside(), 1, i, dim);
                            auto index_i = mapper.subIndex(element, local_index_i, dim);
                            b[index_i] += hcoeff * phiHat[local_index_i] * dy * weight;
                        }
                    }  // end qpoint
                } // end of isNeumann()
                else if(isDirichlet(center)){
                    // Obiđimo sve vrhove na stranici elementa
                    for(int i = 0; i < verticesNumber; i++){
                        auto local_index_i = ref.subEntity(side.indexInInside(), 1, i, dim);
                        auto index_i = mapper.subIndex(element, local_index_i, dim);
                        A[index_i] = 0.0;
                        A[index_i][index_i] = 1.0;
                        auto vertex_i = element.geometry().corner(local_index_i);
                        b[index_i] = exact(vertex_i);
                    }
                } // end of isDirichlet()
            } // end of side.boundary()
        } // end of intersections
    } // end of elements
    // Ispiši matricu radi kontrole.
    //  Dune::writeMatrixToMatlab(A, "matrica");
    std::cout << " Matrica i desna strana su asemblirani.\n";
}
int main(int argc, char *argv[]){
    Dune::MPIHelper::instance(argc, argv);
    const int dim = 3;
    using GridType = Dune::UGGrid<dim>;
    using GridView = typename GridType::LeafGridView;
    // Učitaj mrežu
    GridType *gridptr = Dune::GmshReader<GridType>::read("cube.msh");
    int nref = 0;
    if (argc > 1){
        nref = std::stoi(argv[1]);
    }
    gridptr->globalRefine(nref);
    auto gv = gridptr->leafGridView();
    std::cout << " Učitavanje mreže je gotovo.\n";

    // Linearna algebra
    using Matrix = Dune::BCRSMatrix< Dune::FieldMatrix<double, 1, 1> >;
    using Vector = Dune::BlockVector<Dune::FieldVector<double, 1>    >;

    Matrix A;     // Matrica krutosti
    Vector x, b;  // rješenje i desna strana
    Dune::PkLocalFiniteElement<double, double, dim, 1> fem;

    assemble(gv, fem, A, b);

    // formiranje rješavača
    Dune::MatrixAdapter<Matrix, Vector, Vector> op(A);
    Dune::SeqILU<Matrix, Vector, Vector> ilu1(A, 1, 0.92);
    Dune::BiCGSTABSolver<Vector> bcgs(op, ilu1, 1e-15, 5000, 0);
    Dune::InverseOperatorResult r;

    x.resize(b.N());
    x = 0.0;

    // rješavanje sustava
    bcgs.apply(x, b, r);
    std::cout << " Sustav je riješen.\n";

    if(r.converged){
        std::cout << "Broj iteracija = "        << r.iterations
                  << ", redukcija residuala = " << r.reduction << "\n";
    }
    else{
        std::cout << "Solver nije konvergirao.\n";
    }
    Vector ex;  // vektor točnog rješenja
    vertexdata(gv, ex);

    Dune::VTKWriter<GridView> vtkwriter(gv);
    vtkwriter.addVertexData(x, "u");
    vtkwriter.addVertexData(ex, "exact");
    Vector diff(ex);  // vektor greške
    diff -= x;
    vtkwriter.addVertexData(diff, "diff");
    vtkwriter.write("fem", Dune::VTK::OutputType::ascii);

    std::cout << " Greška aproksimacije (L^infty) = " << diff.infinity_norm()
              << std::endl;

    return 0;
}

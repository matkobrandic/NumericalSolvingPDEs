#ifdef HAVE_CONFIG_H
# include "config.h"     
#endif

#include <dune/common/parallel/mpihelper.hh>

#include <dune/grid/yaspgrid.hh>         
#include <dune/grid/io/file/vtk.hh>      

#include <dune/istl/bvector.hh>
#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/operators.hh>
#include <dune/istl/preconditioners.hh>
#include <dune/istl/solvers.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>

#include <iostream>
#include <cmath>

/*
 * Program za računanje rješenja zadaće 
 *       - div(k grad u ) = f
 * sa Dirichletovim rubnim uvjetom, metodom konačnih volumena.
 * Domena je pravokutna, a koeficijent k je diskontinuiran. U središnjem kvadratu 
 * (bleft,bright)x(bleft,bright) je jednu vrijednost (kInside), a izvan njega 
 * drugu (kOutside). 
 */

// PODACI ZADAĆE -- NE MIJENJATI ----
constexpr int dim = 2;  // dimenzija mreže (mora biti konstantna)

//lijeva i desna granica centralnog bloka
const double bleft = 0.25;
const double bright= 0.75;

// Vrijednosti koeficijenta k u bloku i izvan njega.
const double kInside = 1.0;
const double kOutside = 0.01;

// provjeri da li je varijabla u centralnom bloku (u 1D)
bool check_bounds(double x){
  return (x>bleft) and (x< bright);
}

// provjeri da li je varijabla u centralnom bloku (u dim dimenzija)
template <int dim>
double isInside(Dune::FieldVector<double, dim> const & x){
  bool inside = check_bounds(x[0]);
  for(unsigned int i=1; i<dim; ++i)
     inside = inside and check_bounds(x[i]);
  return inside; 
}

// Koeficijet k diferencijalne jednadžbe 
template <int dim>
double fun_k(Dune::FieldVector<double, dim> const & x){
  if(isInside(x))
     return kInside;
  else
     return kOutside; 
}

/**
  * Dirichletov rubi uvjet.
 */
template <int dim>
double bc(Dune::FieldVector<double, dim> const & x){
    return 0.0;
}

// Desna strana f.
template <int dim>
double fun_f(Dune::FieldVector<double, dim> const & x){
  return 1.0;
}


// "Skalarni" blok-vektor i blok-matrica.
using Vector = Dune::BlockVector<Dune::FieldVector<double,1>>;
using Matrix = Dune::BCRSMatrix<Dune::FieldMatrix<double, 1, 1>>;
// x = vektor rješenja, F = vektor desne strane, coeff = koeficijent k
Vector x,F,coeff,neighbor;
// A je matrica sustava
Matrix A;

//  KRAJ PODATAKA ----------------------------

/**
 * Funkcija koja na danoj mreži formira globani vektor coeff
 * koji sadrži vrijednosti koeficijenta k u centrima elemenata.
 * Koristi se samo za isctravanje u VTK formatu.
 */
template <typename GV>
void calculateK(GV const & gv){
   //auto const & indexSet = gv.indexSet();
   double size = gv.indexSet().size(0);
   coeff.resize(size);
   for (auto const & element : elements(gv)){
      coeff[indexSet.index(element)] = fun_k(element.geometry().center());
   }
}
/** Dimenzioniranje matrice i vektora desne strane.
 *  Računanje profila matrice. 
 */
template <typename GV>
void setMatrixProfile(GV const & gv){
   auto const & elementsNumber = gv.size(0);
   A.setSize(elementsNumber, elementsNumber);
   A.setBuildMode(Matrix::random);
   double size = gv.indexSet().size(0);
   neighbor.resize(size);
   for(auto const & element : elements(gv)){
      //neighbor[gv.indexSet().index(element)]  = 0;
      for(auto const & side : intersections(gv, element)){
         if(side.neighbor())
            neighbor[gv.indexSet().index(element)] += 1;
      }
   }
   for(unsigned int i = 0; i < elementsNumber; ++i){
      A.setrowsize(i, neighbor[i] + 1);
   }
   A.endrowsizes();

   for(unsigned int i = 0; i < elementsNumber; ++i){
      A.addindex(i, i);
      for(unsigned int j = 0; ; ){
         A.addindex(i, j);
      }
   }
   A.addindex(elementsNumber-1, elementsNumber-1);
   A.endindices();
}

/**
 * Asembliranje matrice i vektora desne strane. Dimenzioniranje
 * vektora i profil matrice su obavljeni u setMatrixProfile().
 */
template <typename GV>
void assemble(GV const & gv){
    A = 0.0;
    F = 0.0;

    auto const & is = gv.indexSet();
    for(auto const & E : elements(gv)){
        // VAŠ KOD DOLAZI OVDJE         
     }
}

int main(int argc, char * argv[])
{
    Dune::MPIHelper::instance(argc, argv);

    int Nel = 20;  // pogađa 0.25 i 0.75
    int Nref = 0;  // broj profinjenja
    if(argc > 1) Nref = std::stoi(argv[1]);
    assert(Nref >= 0);

    using GridType = Dune::YaspGrid<dim>;
    using GridView = GridType::LeafGridView;

    Dune::FieldVector<double, dim> L(1.0);             // Duljina stranice
    std::array<int,dim>       s={Nel,Nel};          // broj ćelija po stranici
    // Konstrukcija grida (poziv konstruktora).
    GridType grid(L, s); // serijska mreža
    grid.globalRefine(Nref);
    GridView gv = grid.leafGridView();
      
        // VAŠ KOD DOLAZI OVDJE         
     

   // Ispis rješenja 
    Dune::VTKWriter<GridView> writer(gv);
    writer.addCellData(x,"u");
    // Izračunamo egzaktno rješenje
    calculateK(gv);
    writer.addCellData(coeff, "K");
    writer.write("fvm",Dune::VTK::OutputType::ascii);

    return 0;
}

#include "config.h"
#include <dune/grid/yaspgrid.hh>
#include "dune/common/parallel/mpihelper.hh" 
#include <dune/geometry/quadraturerules.hh>
#include <exception>
#include <iostream>
#include <iomanip>

// parametrizirana funkcija
template<int dim>
double fun(const Dune::FieldVector<double,dim>& x)
{
    return exp(-3.234*(x*x));
}

// Integriraj po elementu sa zadanom točnošću
template<class Entity>
double integrateEntity (const Entity &element, int p)
{
  // dimension of the entity
  const int dim = Entity::dimension;
  // Geometry objekt ima sve informacije o geometriji objekta.
  const auto geometry = element.geometry();
  // Uzmimo sada tip geometrije
  const auto gt = geometry.type();
  // Zatražimo kvadraturnu formulu (obavezno koristiti referencu)
  const auto& rule = Dune::QuadratureRules<double,dim>::rule(gt,p);
  //  provjeri da smo dobili traženi red.
  if (rule.order()<p) std::exit(-1);

  // compute approximate integral
  double result=0;
  for (auto const & qpoint : rule)
    {
      double fval = fun(geometry.global(qpoint.position()));
      double weight = qpoint.weight();
          // | det (grad g) | daje Geometry objekt
      double detjac = geometry.integrationElement(qpoint.position());
      result += fval * weight * detjac;
    }

  return result;
}

template<class GV>
double uniformintegration(GV const& gridView, int p) {
     double value = 0.0;
     for (auto const & element : elements(gridView))
            value += integrateEntity(element, p);
     return value;
}

// Template funkcija. Template parametar je tip mreže.
// Korištenjem template funkcije ne moramo mijenjati funkciju ako u glavnom
// programu promijenimo tip mreže.
// GLV = GridLeafView
template <class GLV> void indices(GLV &leafView) {
  // Uzimamo neke informacije o mreži:
  // Dimenzija mreže
  const int dim = GLV::dimension;
  std::cout << " dim = " << dim << std::endl;
  auto &set = leafView.indexSet();
  std::cout << "Ukupni broj entiteta kodimenzije 0 u mreži = " << set.size(0)
            << std::endl;
  std::cout << "Ukupni broj entiteta kodimenzije 1 u mreži = " << set.size(1)
            << std::endl;
  std::cout << "Ukupni broj entiteta kodimenzije 2 u mreži = " << set.size(2)
            << std::endl;

  int count = 0;
  for (auto const &elem : elements(leafView)) {
    // Geometrijski tip elementa
    auto gt = elem.type();

    // Ispišimo indeks elementa pozivanjem metode index() na elementu.
    std::cout << "Element " << count << "; tip = " << gt
              << ", index = " << set.index(elem) << std::endl;
    // Entity (codim=0) nam može dati broj entiteta veće kodimenzije.
    int n_1 = elem.subEntities(1);
    int n_2 = elem.subEntities(2);
    std::cout << "Broj subentiteta kodimenzije 1 u elementu = " << n_1
              << std::endl;
    for (int i = 0; i < n_1; ++i) {
      std::cout << "Type = " << elem.template subEntity<1>(i).type();
      std::cout << " subIndex = " << set.subIndex(elem, i, 1) << std::endl;
    }

    std::cout << "Broj subentiteta kodimenzije 2 u elementu = " << n_2
              << std::endl;
    for (int i = 0; i < n_2; ++i) {
      std::cout << "Type = " << elem.template subEntity<2>(i).type();
      std::cout << " subIndex = " << set.subIndex(elem, i, 2) << std::endl;
    }

    count++;
  }

  std::cout << "Broj (leaf) elemenata = " << count << std::endl;

  std::cout << std::endl;
}

int main(int argc, char **argv)
{
    Dune::MPIHelper::instance(argc, argv);

    const int dim = 2;
    typedef Dune::YaspGrid<dim> GridType;
    Dune::FieldVector<GridType::ctype, dim> L(2.0);
    std::array<int, dim> N {5,5};
    std::bitset<dim>    periodic(false);
    int overlap = 0;
    GridType grid(L, N, periodic, overlap);
    auto gridView = grid.leafGridView();
    indices(gridView);
    // Integriramo 6 puta i svaki puta profinjujemo mrežu
    double oldvalue = 1E100;
    for (int k = 0; k < 6; k++) {
         double value = uniformintegration(gridView, 2);
         // print result and error estimate
         std::cout << "br elemenata=" << std::setw(8) << std::right << gridView.size(0)
                   << " integral=" << std::scientific << std::setprecision(12)
                   << value << " pocjena greske=" << std::abs(value - oldvalue)
                   << std::endl;

         oldvalue = value;
         grid.globalRefine(1);
    }
    return 0;
}



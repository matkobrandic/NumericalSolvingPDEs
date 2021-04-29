#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <cassert>
#include <cmath>
#include <iomanip>
#include <iostream>

#include <dune/common/fvector.hh>
#include <dune/common/parallel/mpihelper.hh>
#include <dune/geometry/quadraturerules.hh>
#include <dune/grid/yaspgrid.hh>

template <int dim>
using Point = Dune::FieldVector<double, dim>;

template <int dim>
Point<dim> fun(Point<dim> const &x) {
  Point<dim> y(0.0);
  double a = 2.0;
  assert(dim <= 3);

  y[0] = x[0] * x[0] * x[1] * x[1] * x[1] * std::cos(a * (x[0] + x[1]));
  if (dim > 1)
    y[1] = (x[0] * x[0] + x[1] * x[1]) * std::exp(x[0] * x[1]);
  if (dim > 2)
    y[2] = std::sin(a * x[2]) * std::exp(x[1]) * x[2] * x[2] * x[2] * x[2];
  return y;
}

template <int dim> // divergencija
double div(Point<dim> const &x) {
  double a = 2.0;
  double divergencija =
      2 * x[0] * x[1] * x[1] * x[1] * std::cos(a * (x[0] + x[1])) -
      a * x[0] * x[0] * x[1] * x[1] * x[1] * std::sin(a * (x[0] + x[1]));
  if (dim > 1)
    divergencija += 2 * x[1] * std::exp(x[0] * x[1]) +
                    (x[0] * x[0] + x[1] * x[1]) * x[0] * std::exp(x[0] * x[1]);
  if (dim > 2)
    divergencija +=
        a * std::cos(a * x[2]) * std::exp(x[1]) * x[2] * x[2] * x[2] * x[2] +
        4 * std::sin(a * x[2]) * std::exp(x[1]) * x[2] * x[2] * x[2];

  return divergencija;
}

/* Integracija po rubu domene vektorskog polja f.n.
 * p = red točnosti integracijske formule
 */
template <typename Grid>
double bdry_integration(Grid &grid, int p) {

  const int dim = Grid::dimension;
  auto gridView = grid.leafGridView();
  double integral = 0.0;
  // Petlja po svim elementima
  for (auto const &element : elements(gridView)) {
    double elem_integral = 0.0;
    // petlja po svim stranicama elementa
    for (auto const &side : intersections(gridView, element)) {

      if (side.boundary()) // Jesmo li na granici domene?
      {
        const auto sidegeo = side.geometry();
        auto outerNormal = side.centerUnitOuterNormal();
        // Zatražimo kvadraturnu formulu na stranici (dimenzija dim -1)
        const auto &rule =
            Dune::QuadratureRules<double, dim - 1>::rule(sidegeo.type(), p);
        if (rule.order() < p)
             std::cerr << "Integracijska formula reda " << p << " nije dostupna.\n";
        double result = 0.0;
        // Petlja po svim integracijskim točkama
        for (auto const &qpoint : rule) {
          auto fval = fun(sidegeo.global(qpoint.position()));
          double weight = qpoint.weight();
          // | det (grad g) | daje Geometry objekt
          double detjac = sidegeo.integrationElement(qpoint.position());
          result += (fval * outerNormal) * weight * detjac;
        }

        elem_integral += result;
      }
    } // kraj petlje po svim stranicama
    integral += elem_integral;
  } // kraj petlje po svim elementima

  return integral;
}


int main(int argc, char *argv[]) {
  Dune::MPIHelper::instance(argc, argv);

  const int dim = 2;
  typedef Dune::YaspGrid<dim> GridType;
  Dune::FieldVector<GridType::ctype, dim> L(3.0);
  std::array<int, dim> N{5, 5};
  std::bitset<dim> periodic(false);
  int overlap = 0;
  GridType grid(L, N, periodic, overlap);

  for (int k = 1; k < 5; ++k) {
    grid.globalRefine(1);
    double bdry_integral = bdry_integration(grid, 5);
    std::cout << "elements=" << std::setw(8) << std::right << grid.size(0)
              << " bdry_integral = " << std::scientific << std::setprecision(12)
              << bdry_integral << "\n";
  }

  return 0;
}

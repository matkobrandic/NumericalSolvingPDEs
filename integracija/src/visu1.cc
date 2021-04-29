#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <iostream>

#include <dune/common/fvector.hh>
#include <dune/grid/common/mcmgmapper.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/grid/yaspgrid.hh>

#include "elementdata.hh"

// Funkcija koju prikazujemo
template <int dim>
//double f(Dune::FieldVector<double, 3> const &x) { return 2 * (x * x); }
double f(Dune::FieldVector<double, dim> const &x) { return 2 * (x * x); }

int main(int argc, char **argv) {
  Dune::MPIHelper::instance(argc, argv);

  // Broj profinjenja mreže
  int refSteps = 2;
  // Argument komandne linije može biti broj profinjenja mreže
  if (argc > 1) {
    refSteps = std::stoi(argv[1]); // stoi -> konvertira string u int
  }

  //const int dim = 3;
  const int dim = 2;
  typedef Dune::YaspGrid<dim> GridType;
  typedef GridType::LeafGridView GridView;
  Dune::FieldVector<GridType::ctype, dim> L(1.0);
  //std::array<int, dim> s = {1, 1, 1};
  std::array<int, dim> s = {1, 1};
  //std::bitset<dim> periodic(false); //jer default vrijednosti će biti sasvim dobre
  int overlap = 0;
  //GridType grid(L, s, periodic, overlap);
  GridType grid(L, s);
  GridView leafView = grid.leafGridView();

  std::string filename = "data";
  // Instanciranje funkcije
  //  Exp<double> f;

  // Profini mrežu
  grid.globalRefine(refSteps);
  // Vektor koeficijenata
  std::vector<double> coeff;

  // Izračunaj coeff
  elementdata(leafView, f<2>, coeff);

  // Ispisivanje podataka u VTK datoteku.
  Dune::VTKWriter<GridView> vtkwriter(leafView);
  // Imamo metode addCellData i addVertexData za dodavanje podataka vezanih uz
  // elemente odnosno vrhove. "data" je ime podatka.
  vtkwriter.addCellData(coeff, "el_data");
  vtkwriter.write("element_" + filename, Dune::VTK::OutputType::ascii);

  return 0;
}

/** \file

    \brief Parabolička jednadžba diskretizirana konformnom metodom KE u prostoru
          i implicitnom Eulerovom metodom u vremenu.
    */
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#include <array>
#include <bitset>
#include <cmath>
#include <iostream>
#include <map>
#include <string>
#include <vector>

#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>
#include <dune/common/parallel/mpihelper.hh>
#include <dune/grid/io/file/gmshreader.hh>
#include <dune/grid/io/file/vtk/subsamplingvtkwriter.hh>
#include <dune/grid/yaspgrid.hh>
#if HAVE_UG
#include <dune/grid/uggrid.hh>
#endif
#if HAVE_ALUGRID
#include <dune/grid/alugrid.hh>
#include <dune/grid/io/file/dgfparser/dgfalu.hh>
#include <dune/grid/io/file/dgfparser/dgfparser.hh>
#endif

#include "driver.hh"
//#include <fenv.h>

int main(int argc, char **argv){
	// Inicijaliziraj Mpi
	Dune::MPIHelper &helper = Dune::MPIHelper::instance(argc, argv);

	if (argc != 4) {
		if (helper.rank() == 0)
			std::cout << "usage: ./parabolic <refine level> <init dt> <tend>" << std::endl;
		return 1;
	}
	//  feenableexcept(FE_INVALID | FE_OVERFLOW);
	int    level   = std::stoi(argv[1]);  // broj profinjenja mreže
	double dt      = std::stod(argv[2]);  // vremenski korak
	double tend    = std::stod(argv[3]);  // krajnje vrijeme simulacije

	// Konstrukcija mreže
	constexpr int dim = 2;
	Dune::FieldVector<double, dim> ll{-100.0, -100.0};
	Dune::FieldVector<double, dim> ur{100.0, 100.0};
	std::array<int, dim> N{10, 10};
	using coo = Dune::EquidistantOffsetCoordinates<double, dim>;
	Dune::YaspGrid<dim, coo> grid(ll, ur, N);
	grid.globalRefine(level);
	const auto &gv = grid.leafGridView();
	driver(gv, dt, tend);
}

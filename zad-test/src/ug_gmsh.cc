#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

#include <iostream>

#include <dune/grid/uggrid.hh>             // Koristimo UGGrid
#include <dune/grid/io/file/gmshreader.hh> // GmshReader klasa
#include <dune/grid/common/gridinfo.hh>
#include <dune/common/parallel/mpihelper.hh>
#include <dune/grid/io/file/vtk.hh>

int main(int argc, char** argv)
{
	if (argc < 2) {
		std::cerr << "Usage: " << argv[0] << " ime_grid_datoteke.msh"
				<< std::endl;
		std::exit(1);
	}
	Dune::MPIHelper& helper = Dune::MPIHelper::instance(argc, argv);

    using GridType = Dune::UGGrid<2>;
    using GridView = GridType::LeafGridView;

    bool verbosity = true;
    bool insertBoundarySegments = false;  // Bez toga Dune::GmshReader zna podbaciti (barem u 3D)
    
    GridType* pgrid = Dune::GmshReader<GridType>::read(argv[1], verbosity, insertBoundarySegments);
    // loadBalance() distribuira mrežu po procesorima
    pgrid->loadBalance();

    Dune::gridinfo(*pgrid);
    Dune::VTKWriter<GridView> writer(pgrid->leafGridView());
    writer.write("ug_gmsh");

	return 0;
}





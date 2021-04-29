#ifdef HAVE_CONFIG_H
#include "config.h"     
#endif
#include <iostream>
#include <string>

#include <dune/grid/common/gridinfo.hh> 
// Za ALUGrid
#include <dune/alugrid/grid.hh>
// Za čitanje ALUGrida iz DGF datoteke trebamo dvije datoteke          
#include <dune/grid/io/file/dgfparser/gridptr.hh>
#include <dune/alugrid/dgf.hh>
#include <dune/grid/io/file/vtk.hh>


int main(int argc, char** argv) {
    Dune::MPIHelper::instance(argc, argv);
  
    // trebamo barem jedan argument -- ime datoteke
    if(argc < 2) {
	    std::cerr << "Usage: " << argv[0] << " grid_name.dgf [no_refinement]" << std::endl;
	    std::exit(1);
    }	 
    std::string ime_dat(argv[1]); // ime datoteke je dano u prvom argumentu
    // Pretvori string u int ili izbaci izuzetak.
    int no_r = 0; // nema profinjenja ako drugi argument nije zadan
    if(argc == 3)
        no_r = std::stoi(argv[2]);  // stoi (= string to int)
    
    const int dim = 2;

    using GridType = Dune::ALUGrid<dim, dim, Dune::cube, Dune::nonconforming>;
    using GridView = GridType::LeafGridView;

    Dune::GridPtr<GridType> gridptr(ime_dat);
    GridType & grid = *gridptr;
    // podijeli svaki četverokut na 4  dijela no_r puta
    grid.globalRefine(no_r);

    // print some information about the grid
    Dune::gridinfo(grid);
    Dune::VTKWriter<GridView> writer(grid.leafGridView());
    writer.write("alu_dgf");

    return 0;
}

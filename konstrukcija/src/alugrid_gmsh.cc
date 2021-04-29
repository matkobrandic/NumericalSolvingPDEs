#ifdef HAVE_CONFIG_H
# include "config.h"     // Datoteka koju kreira configure skripta i koja služi
#endif                   // adaptaciji koda na okolinu u kojoj se kompilira. 

#include <iostream>     // Ulazno-izlazna biblioteka 
#include <cstdlib>     

#include "dune/common/parallel/mpihelper.hh"  // Inicijalizacija MPI sustava (za paralelno izvršavanje programa)

#include <dune/alugrid/grid.hh>           // koristimo ALUGrid
#include <dune/grid/io/file/gmshreader.hh> // GmshReader klasa
#include <dune/grid/common/gridinfo.hh>    // Informacije o mreži. 
#include <dune/grid/io/file/vtk.hh>

// glavni program 
// Ime datoteke tipa .msh očekujemo kao prvi argument komandne linije.
int main(int argc, char** argv)
{
    Dune::MPIHelper::instance(argc, argv);
  
    if(argc < 3){
       std::cerr << "Usage: " << argv[0] << " ime_grid_datoteke.msh br_profinjenja" << std::endl;
       std::exit(1);
    }

    using GridType =  Dune::ALUGrid<2,2,Dune::simplex,Dune::conforming>;
    using GridView = GridType::LeafGridView;
    // statička metoda  Dune::GmshReader<GridType>::read vraća pokazivač na kreirani grid.
    GridType* pgrid = Dune::GmshReader<GridType>::read(argv[1]);
    int no_r = std::stoi(argv[2]);
    pgrid->globalRefine(no_r);     // profini mrežu
    // Ispisuje neke informacije o mreži
    Dune::gridinfo(*pgrid);
    Dune::VTKWriter<GridView> writer(pgrid->leafGridView());
    writer.write("alu_gmsh");
    return 0;
}

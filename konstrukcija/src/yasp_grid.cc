#ifdef HAVE_CONFIG_H
# include "config.h"     // Datoteka koja služi
#endif                   // adaptaciji koda na okolinu u kojoj se kompilira. 
                         // Kada se ispusti znaju se javljati "čudne greške". 

#include <iostream>     // Ulazno-izlazna biblioteka

#include <dune/grid/yaspgrid.hh>         // Koristimo YaspGrid
#include <dune/grid/common/gridinfo.hh>  // za funkciju gridinfo()
#include <dune/grid/io/file/vtk.hh>      // za VTKWriter

// glavni program 
int main(int argc, char** argv)
{
//  Čak i u sekvencijalnom programu moramo pozvati MPIHelper.   
    Dune::MPIHelper::instance(argc, argv);

    const int dim = 3;  // dimenzija mreže (mora biti konstantna)
    // U većim programima je korisno uvesti kraća imena za tipove pomoću using ili typedef naredbe.
    using GridType = Dune::YaspGrid<dim>;
    using GridView = GridType::LeafGridView;

    // FieldVector je vektor konstantne duljine.
    Dune::FieldVector<double, dim> L(1.0);             // Duljina stranice
    L[1]= 2.0;
    std::array<int,dim>            s={16,26,16};          // broj ćelija po stranici
    std::bitset<dim>               periodic(false);    // periodičnost u svim smjerovima
    int overlap = 0;                                   // preklapanje za paralelni grid
    // Konstrukcija grida (poziv konstruktora).
    GridType grid(L, s, periodic, overlap); // serijska mreža
    // printa neke informacije o mreži
    Dune::gridinfo(grid);

    Dune::VTKWriter<GridView> writer(grid.leafGridView()); // Ispis mreže u
    writer.write("yasp");                                  // datoteku yasp.vtu

    return 0;
}

#ifdef HAVE_CONFIG_H
# include "config.h"     
#endif
#include <iostream>
#include <vector>

#include <dune/alugrid/grid.hh>  
#include <dune/grid/io/file/dgfparser/gridptr.hh>
#include <dune/alugrid/dgf.hh>
// Za ispis mreže u VTK formatu
#include <dune/grid/io/file/vtk/vtkwriter.hh>


int main(int argc, char** argv)
{
    Dune::MPIHelper::instance(argc, argv);
    
    const int dim = 3;
    using GridType = Dune::ALUGrid<dim,dim,Dune::simplex,Dune::conforming>;
    using GridView = GridType::LeafGridView ;  // Treba nam za  Dune::VTKWriter

    if(argc < 2){
      std::cout << "Need  one argument - a msh file name!" << std::endl;
      std::exit(-1);
    } 

    // Kreiramo dvodimenzionalnu mrežu čitanjem gmsh datoteke
    //GridType * p_grid = Dune::GmshReader<GridType>::read(argv[1]);
    Dune::GridPtr<GridType> p_grid(argv[1]);
    // Uzmi GridView.
    auto gridView = p_grid->leafGridView();

    // Dimenzija mreže
    std::cout << " dim = " << dim << std::endl;
   
    // Tip globalne koordinate
    using GlobalCoordinate = GridView::template Codim<0>::Geometry::GlobalCoordinate;

    // Iteriramo po svim elementima mreže i :
    // 1) brojimo sve elemente
    // 2) računamo volumen domene 
    // 3) ispisujemo sve vrhove elemenata 
    // 4) ispisujemo tip elementa
    // 5) ispisujemo centar elementa
    int count = 0;
    double volumen = 0.0;
    for(auto const & element : elements(gridView))
    {
        // Geometrijski tip elementa (tip Dune::GeometryType)
        auto gt = element.type();
        auto geom = element.geometry();
        // Broj vrhova elementa
        int n_v = geom.corners();

        // Uzmi koordinate svih vrhova elemeta
        std::vector<GlobalCoordinate> coo_v(n_v);
        //std::vector<decltype(geom.corner(0))> coo_v(n_v); // drugi način detekcije tipa koordinate
        for(unsigned int i=0; i <n_v; ++i) coo_v[i] = geom.corner(i);

        // volumen elementa
        double vol = geom.volume();
        // Koordinate centra elementa
        auto centar = geom.center();

        std::cout <<"Element " << count << "; tip = " <<  gt << ", volumen = " << vol
                  << "\n Koordinate vrhova: ";
        for(unsigned int i=0; i <n_v; ++i) std::cout << coo_v[i] <<",";
        std::cout <<"\n Centar elementa: " << centar << std::endl;

        count++; 
        volumen += vol;
    } 

    std::cout << "Broj (leaf) elemenata = " << count  << std::endl;
    std::cout << "Volumen domene  = " << volumen  << std::endl;

    std::cout << std::endl;

    Dune::VTKWriter<GridView> vtkwriter(gridView);
    vtkwriter.write("output");

    return 0;
}

#ifdef HAVE_CONFIG_H
# include "config.h"     
#endif
#include <iostream>
#include <vector>
#include <cassert>
#include <cmath>
#include <limits>

#include "dune/common/parallel/mpihelper.hh"
#include <dune/common/exceptions.hh> 

#include <dune/grid/uggrid.hh>  
#include <dune/grid/common/gridinfo.hh> 
#include <dune/grid/io/file/gmshreader.hh>

#include <dune/grid/io/file/vtk/vtkwriter.hh>

// Izračunaj kut između p1-p2 i p3-p2
template <typename Point>
double kut(Point const & p1, Point const & p2, Point const & p3)
{
    auto vec_one = p1 - p2;
    auto vec_two = p3 - p2;

    return acos(vec_one.dot(vec_two) / (vec_one.two_norm() * vec_two.two_norm()));
}

int main(int argc, char** argv)
{
    const int dim = 2;
    typedef Dune::UGGrid<dim> GridType;
    GridType * p_grid =nullptr;

    Dune::MPIHelper& helper = Dune::MPIHelper::instance(argc, argv);

    if(argc < 2){
      std::cout << "Need msh file name!" << std::endl;
      std::exit(-1);
    } 

    // Kreiramo dvodimenzionalnu mrežu čitanjem gmsh datoteke
    p_grid = Dune::GmshReader<GridType>::read(argv[1]);

    // Ako je mreža paralelna rasporedi je ravnomjerno po procesorima
    p_grid->loadBalance();
    if(helper.rank() == 0) Dune::gridinfo(*p_grid);
   
    // Uzmi GridView 
    typedef typename GridType::LeafGridView LeafGridView;
    LeafGridView gridView = p_grid->leafGridView(); 
    
    double min = std::numeric_limits<double>::max();  // najveći double
    double max = std::numeric_limits<double>::lowest(); // najveći negativni double
    int count = 0;
    double loc_min;
    double loc_max;
    for(auto const & element : elements(gridView))
    {
        //getting the points of an element
		int elem_number = element.geometry().corners();
        std::vector<Dune::FieldVector<double,dim> > p(elem_number);

        for(unsigned int i = 0; i < elem_number; ++i)
			p[i] = element.geometry().corner(i);

        //for(int i = 1; i < elem_number-1; ++i)
            //getting the angles
        /*
            double angle[elem_number];
            //error, angle nije definiran (???)
            angle[0] = kut(p[elem_number-1], p[0], p[1]);
            loc_min = angle[0];
            loc_max = angle[0];
            double angle_current;
            for(int i = 1; i < elem_number-1; ++i){
                angle_current = kut(p[i-1], p[i], p[i+1]);
                angle[i] = angle_current;
                if(angle_current < loc_min)
                    loc_min = angle[i];
                if(angle_current < loc_max)
                    loc_max = angle[i];
            }
            angle[elem_number-1] = kut(p[elem_number-2], p[elem_number-1], p[0]);
            if(angle[elem_number-1] < loc_min)
                loc_min = angle[elem_number-1];
            if(angle[elem_number-1] < loc_max)
                loc_max = angle[elem_number-1];
        */
        loc_min = kut(p[elem_number-1], p[0], p[1]);
        loc_max = loc_max;
        double angle_current;
        for(int i = 1; i < elem_number-1; ++i){
            angle_current = kut(p[i-1], p[i], p[i+1]);
            if(angle_current < loc_min)
                loc_min = angle_current;
            if(angle_current > loc_max)
                loc_max = angle_current;
        }
        angle_current = kut(p[elem_number-2], p[elem_number-1], p[0]);
        if(angle_current < loc_min)
            loc_min = angle_current;
        if(angle_current > loc_max)
            loc_max = angle_current;
        //ispis najvećeg i najmanjeg elementa za svaki element mreže
        std::cout <<"Element " << count << " min kut = " <<  180*loc_min/M_PI
                << ", max kut = " << 180*loc_max/M_PI << "\n";

        min = std::min(loc_min, min);
        max = std::max(loc_max, max);
        count++;
    }

    std::cout << "Broj (leaf) elemenata = " << count  << std::endl;
    std::cout << "Minimalni kut  = " << 180*min/M_PI 
              << ", maksimalni kut = " << 180*max/M_PI  << std::endl;


    // Ispis mreže u VTK formatu (u datoteku poluvijenac.vtu)
    Dune::VTKWriter<LeafGridView> vtkwriter(gridView);
    vtkwriter.write("poluvijenac");

    delete p_grid;
 
    return 0;
}

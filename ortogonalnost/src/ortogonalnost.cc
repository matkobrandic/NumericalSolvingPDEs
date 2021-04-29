#ifdef HAVE_CONFIG_H
# include "config.h"
#endif
#include <iostream>
#include <cmath>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/exceptions.hh>
#include <dune/grid/uggrid.hh>
#include <dune/grid/io/file/gmshreader.hh>
#include <dune/grid/io/file/vtk.hh>
#include <dune/grid/common/gridinfo.hh>

template<typename point>
double angle(point p1, point p2, point p3,point p4){
    point first = p1 - p2;
    point second = p3 - p4;
    double angle = acos(first.dot(second) / (first.two_norm() * second.two_norm()));
    if (angle > M_PI/2){
        angle = M_PI - angle;
    }
    return M_PI/2 - angle;
}
// GV = GridView type
template <typename GV>
double ortogonalnost(GV const & gridview){
    double delta = 0;
    for(const auto& element : elements(gridview)){
        for(const auto& side : intersections(gridview, element)){
            //auto geom = side.geometry();
            auto p1 = element.geometry().center();
            decltype(p1) p2; 
            if (side.neighbor()){
                p2 = side.outside().geometry().center();
            }
            else{
                p2 = side.geometry().center();
            }
            auto p3 = side.geometry().corner(0);
            auto p4 = side.geometry().corner(1);
            delta = std::max(delta,angle(p1,p2,p3,p4));
        }
    }
    return delta * 180 / M_PI;
}
int main(int argc, char** argv){
    //initialize
    const int dim = 2;
    typedef Dune::UGGrid<dim> GridType;
    typedef GridType::LeafGridView GridView;
    Dune::MPIHelper& helper = Dune::MPIHelper::instance(argc, argv);

    if(argc < 2){
        std::cout << "Need .msh file!" << std::endl;
        std::exit(-1);
    }

    GridType * grid = Dune::GmshReader<GridType>::read(argv[1]);
    grid->loadBalance();
    if(helper.rank() == 0){
        Dune::gridinfo(*grid);
    }
    auto deviation = ortogonalnost(grid->leafGridView());
    std::cout << "\nDeviation from orthogonality: " << deviation << " degrees.\n";

    Dune::VTKWriter<GridView> vtkwriter(grid->leafGridView());
    vtkwriter.write("poluvijenac");
    
    delete grid;
    return 0;
}


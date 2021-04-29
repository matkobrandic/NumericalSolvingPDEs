#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

#include <iostream>
#include <dune/common/parallel/mpihelper.hh> // An initializer of MPI
#include <dune/common/exceptions.hh>         // We use exceptions
#include <dune/grid/uggrid.hh>               // Koristimo UGGrid
#include <dune/grid/io/file/gmshreader.hh>   // GmshReader klasa
#include <dune/grid/io/file/vtk.hh>
#include <dune/geometry/quadraturerules.hh>

#include <string>
#include <array>

template<int dim>
double fun_xx(const Dune::FieldVector<double,dim>& x){
    return x[0] * x[0];
}
template<int dim>
double fun_y(const Dune::FieldVector<double,dim>& x){
    return x[1];
}
template<int dim>
double fun_x(const Dune::FieldVector<double,dim>& x){
    return x[0];
}
template<int dim>
double fun_1(const Dune::FieldVector<double,dim>& x){
    return 1;
}
template<class Entity, class fun>
double integrateEntity(const Entity &element,int p,fun templateFunction){
    //initialize
    const int dim = Entity::dimension;
    const auto geometry = element.geometry();
    const auto geometryType = geometry.type();
    const auto& rule = Dune::QuadratureRules<double,dim>::rule(geometryType,p);
    //integrating
    double result = 0;
    for(const auto& point : rule){
        double temp = templateFunction(geometry.global(point.position()));
        double weight = point.weight();
        double gradient = geometry.integrationElement(point.position());
        result += temp * weight * gradient;
    }
    return result;
}
template <typename GV>
std::array<double,3> tromost(GV const & gridView){
    //initialize
    std::array<double,3> value = {0.0,0.0,0.0};
    double area = 0;
    //computing
    for(const auto& element : elements(gridView)){
        const auto geom = element.geometry();
        value[0] += integrateEntity(element,2,fun_x<2>);
        value[1] += integrateEntity(element,2,fun_y<2>);
        value[2] += integrateEntity(element,2,fun_xx<2>);
        area += integrateEntity(element,2,fun_1<2>);
    }
    value[0] = value[0]/area;
    value[1] = value[1]/area;
    return value;
}
int main(int argc, char** argv){
    //Initialize
    Dune::MPIHelper& helper = Dune::MPIHelper::instance(argc, argv);
    const int dim = 2;
    using GridType = Dune::UGGrid<dim>;
    using GridView = GridType::LeafGridView;
    bool verbosity = true;
    bool insertBoundarySegments = false;  // Bez toga Dune::GmshReader zna podbaciti (barem u 3D)
    // mashgrid filename
    std::string file("pinch-2D-simplex.msh"); 
    // Read mashgrid from msh file
    GridType* pgrid = Dune::GmshReader<GridType>::read(file, verbosity, insertBoundarySegments);
    GridView gv = pgrid->leafGridView();
    //In case of parallel computing -> distribute grid 
    pgrid->loadBalance();
    // Computing
    auto solution = tromost(gv);
    std::cout << "Težište = (" << solution[0] << "," << solution[1] 
              << "). Tromost = " << solution[2] << std::endl;
    
    Dune::VTKWriter<GridView> writer(gv);
    writer.write("pinch-2D-simplex");

    return 0;
}

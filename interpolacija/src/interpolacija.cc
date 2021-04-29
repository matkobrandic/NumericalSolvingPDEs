#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

#include <cassert>
#include <cmath>
#include <numeric>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/grid/yaspgrid.hh>
#include <dune/grid/io/file/vtk.hh>
#include <dune/localfunctions/lagrange/qk.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>

// Nova verzija Dune-a zahtjeva da eksportiramo DomainType i RangeType.
// Bez toga funkcija interpolate ne radi.
template <int dim>
struct ScalarFunctionTraits{
    typedef typename Dune::FieldVector<double, dim> DomainType;
    typedef typename Dune::FieldVector<double, 1>   RangeType;
};
// Funkcija koju treba interpolirati.
template <typename Element, int dim>
struct Function{
    using Traits = ScalarFunctionTraits<dim>;
    Element const & element;
    Function(Element const & el) : element(el) {}

    void evaluate(typename Traits::DomainType const & x_local,
                  typename Traits::RangeType  & y) const {
        auto x_global = element.geometry().global(x_local);
        evaluate_global(x_global, y);
    }

    void evaluate_global(typename Traits::DomainType const & x_global,
                         typename Traits::RangeType & y) const {
        double x =x_global[0];
	for(unsigned int i=1; i<dim; ++i) // računa x_0*x_1*...*x_{d-1}
        x *= x_global[i];
        y = std::sin(2*M_PI*x);
    }
};
// Računanje interpolacije funkcije. Metoda je neovisna o stupnju
// polinoma koji koji je iskorišten.
int main(int argc, char * argv[]){
    Dune::MPIHelper::instance(argc, argv);
    // Grid
    constexpr int dim = 2;
    // Konstrukcija YaspGrid mreže
    int Nel = 5;
    using GridType = Dune::YaspGrid<dim>;
    using GridView = GridType::LeafGridView;

    Dune::FieldVector<double, dim> L(1.0);
    std::array<int,dim> s = {Nel,Nel};
    GridType grid(L, s);
    // Profinjenje
    grid.globalRefine(4);
    const auto & gridview = grid.leafGridView();
    int verticesNumber = gridview.indexSet().size(2);

    // Konačni element tipa Qk. Eksperimentirati sa k=1 i k=2.
    constexpr int k = 2;
    using FEM = Dune::QkLocalFiniteElement<double, double, dim, k>;
    using DomainType = FEM::Traits::LocalBasisType::Traits::DomainType;
    using RangeType = FEM::Traits::LocalBasisType::Traits::RangeType;
    using Element = GridView::template Codim<0>::Entity;
    FEM fem;

    auto const & basis = fem.localBasis();
    auto const & interpolation = fem.localInterpolation();
    std::vector<double> interpolated(verticesNumber);
    std::vector<double> exact(verticesNumber);

    for(auto element : elements(gridview)){
        auto macroElement = element.father().father().father().father();
        Function<Element, dim> f(macroElement);
        std::vector<double> coeffitients;
        interpolation.interpolate(f, coeffitients);
        for(int i = 0; i < 4; i++){
            auto macroIndex = gridview.indexSet().subIndex(element, i, 2);
            auto x = element.geometry().corner(i);
            RangeType y;
            f.evaluate_global(x, y);
            exact[macroIndex] = y[0];
            std::vector<RangeType> phi;
            basis.evaluateFunction(x,phi);
            interpolated[macroIndex] = 0.0;
            for(int j = 0; j < phi.size(); j++){
                interpolated[macroIndex] += coeffitients[j] * phi[j];
            }
        }
    }

    Dune::VTKWriter<GridView> vtkwriter(gridview);
    vtkwriter.addVertexData(interpolated, "Interpolated");
    vtkwriter.addVertexData(exact, "Exact");
    vtkwriter.write("interpolacija");

    return 0;
}

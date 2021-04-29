#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

#include <array>
#include <bitset>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/grid/yaspgrid.hh>
#include<dune/common/parametertreeparser.hh>

#include "driver.hh"
#include "driver_simple_iter.hh"

int main(int argc, char** argv){
	Dune::MPIHelper::instance(argc, argv);

	// otvaranje datoteke s parametrima - inicijalizacija
	Dune::ParameterTree ptree;
	Dune::ParameterTreeParser ptreeparser;
	ptreeparser.readINITree("nlin.ini",ptree); // datoteka se zove nlin.ini
	ptreeparser.readOptions(argc,argv,ptree);

	// read ini file
	const double  H = ptree.get<double>("H");
	const int level = ptree.get<int>("refinement");
	constexpr int dim = 2;

	// sekvencijalna verzija -- kreiraj Grid
	Dune::FieldVector<double,dim> L(1.0);
	std::array<int,dim>           N{10,10};
	Dune::YaspGrid<dim> grid(L,N);

	grid.globalRefine(level);
	const auto& gv=grid.leafGridView();

	std::cout << "Rješavanje Newtonovom metodom. ===========\n";
	driver(gv, H);
	std::cout <<  "Rješavanje sekvencijalnom metodom. ==========\n";
	driver_simple_iter(gv, H);

	return 0;
}

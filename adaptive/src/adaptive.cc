// Adaptivna metoda konačnih elemenata. Primjer je uzet iz
// dune-pdelab-tutorials/tutorial05.
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <cmath>
#include <iostream>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/parametertreeparser.hh>

#include <dune/grid/uggrid.hh>

#include <dune/alugrid/dgf.hh>
#include <dune/alugrid/grid.hh>
#include <dune/grid/io/file/dgfparser/dgfparser.hh>
#include <dune/grid/io/file/dgfparser/gridptr.hh>
#include <dune/grid/io/file/dgfparser/dgfug.hh>

#include "driver.hh"

int main(int argc, char **argv){
	Dune::MPIHelper::instance(argc, argv);
	const int dim = 2;

	// čitanje domene i parametara
	Dune::ParameterTree ptree;
	Dune::ParameterTreeParser ptreeparser;
	ptreeparser.readINITree("adaptive.ini", ptree);
	ptreeparser.readOptions(argc, argv, ptree);
	std::string filename = "src_dir/unitcube.dgf"; //Učitavamo datoteku za domenu (0,1)x(0,1)
	int subsampling = ptree.get<int>("output.subsampling", 1);
	std::string output = ptree.get<std::string>("output.filename", "output");
	double  tol = ptree.get<double>("fem.tol", 1E-3);
	double alpha = ptree.get<double>("fem.alpha", 0.5);
	double beta = ptree.get<double>("fem.beta", 0.1);
	int steps = ptree.get<int>("fem.steps", 5);

	//error throws
	if(alpha > 1.0 or alpha <0)
		throw std::runtime_error("alpha out of bounds!");
	if(beta > 1.0 or beta <0)
		throw std::runtime_error("beta out of bounds!");
	if(alpha+beta > 1.0)
		throw std::runtime_error("alpha + beta > 1!");
	// Primjer korištenja UG Grida.
	using Grid = Dune::UGGrid<dim>;
	Dune::GridPtr<Grid> gridp(filename);
	driver<Grid>(*gridp, subsampling,  steps, alpha, beta, tol, output);
	return 0;
}
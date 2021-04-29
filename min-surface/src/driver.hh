#ifndef DRIVER_HH
#define DRIVER_HH

#include <dune/istl/bvector.hh>
#include <dune/istl/io.hh>
#include <dune/istl/operators.hh>
#include <dune/istl/preconditioners.hh>
#include <dune/istl/solvers.hh>
#include <dune/istl/superlu.hh>

#include <dune/pdelab/backend/istl.hh>
#include <dune/pdelab/common/function.hh>
#include <dune/pdelab/common/vtkexport.hh>
#include <dune/pdelab/constraints/conforming.hh>
#include <dune/pdelab/finiteelementmap/qkfem.hh>
#include <dune/pdelab/gridfunctionspace/genericdatahandle.hh>
#include <dune/pdelab/gridfunctionspace/gridfunctionspace.hh>
#include <dune/pdelab/gridfunctionspace/gridfunctionspaceutilities.hh>
#include <dune/pdelab/gridfunctionspace/interpolate.hh>
#include <dune/pdelab/gridoperator/gridoperator.hh>
//#include <dune/pdelab/stationary/linearproblem.hh>
#include<dune/pdelab/newton/newton.hh>

#include <dune/grid/io/file/gmshreader.hh>
#include <dune/grid/io/file/vtk/subsamplingvtkwriter.hh>

#include "bctype.hh"
#include "operator.hh"

	
template<typename GV>
void driver (const GV& gv, double H){
	using namespace Dune::PDELab;  // Da skratimo imena

	// Prostor konačnih elemenata (grid function space) i GridOperator.
	const int dim = GV::dimension;
	const int k = 1;
	typedef QkLocalFiniteElementMap<GV,double,double,k>               FEM;
	typedef ConformingDirichletConstraints                            CONSTRAINTS;
	typedef ISTL::VectorBackend<>                                     VBE;
	typedef GridFunctionSpace<GV,FEM,CONSTRAINTS,VBE>                 GFS;
	typedef typename GFS::template ConstraintsContainer<double>::Type CC;
	typedef MinimalSurfaceLOP<DirichletBdry, FEM>                     LOP;
	typedef ISTL::BCRSMatrixBackend<>                                 MBE;
	typedef GridOperator<GFS,GFS,LOP,MBE,double,double,double,CC,CC>  GO;

	FEM fem(gv);
	GFS gfs(gv,fem);
	CC cc;
	DirichletBdry bctype;
	constraints(bctype, gfs, cc); // asembliranje ograničenja Dirichletovog tipa
	std::cout << "constrained dofs=" << cc.size() << " of " << gfs.globalSize() << std::endl;
	LOP lop(H, bctype, fem);
	MBE mbe(9);  // traži prosječan broj ne-nul elemenata u redu (=9)
	GO go(gfs,cc,gfs,cc,lop,mbe);

	//  Konstrukcija rješavača.
	typedef typename GO::Traits::Domain            U;  // U = Dune::PDELab::Backend::Vector<GFS,RF>;
	typedef BCExtension<GV,double>                 G;
	typedef ISTLBackend_SEQ_BCGS_SSOR              LS;

	U u(gfs,0.0); // U = Dune::PDELab::ISTL::BlockVector
	G g(gv);
	interpolate(g, gfs, u);
	LS ls(5000,true);  // max 5000 iteracija, verbosity = true

	Dune::PDELab::Newton<GO,LS,U> newton(go,u,ls);
	newton.setReassembleThreshold(0.0);
	newton.setVerbosityLevel(3);
	newton.setReduction(1e-10);
	newton.setMinLinearReduction(1e-4);
	newton.setMaxIterations(25);
	newton.setLineSearchMaxIterations(30);
	newton.apply();

	// Grafički izlaz (VTK)
	typedef DiscreteGridFunction<GFS,U> DGF;

	DGF udgf(gfs,u);
	Dune::VTKWriter<GV> vtkwriter(gv,Dune::VTK::conforming);
	vtkwriter.addVertexData(std::make_shared<VTKGridFunctionAdapter<DGF>>(udgf,"elevation"));
	vtkwriter.write("surface", Dune::VTK::ascii);

}

#endif

#ifndef DRIVER_SIMPLE_ITER_HH
#define DRIVER_SIMPLE_ITER_HH

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
#include <dune/pdelab/stationary/linearproblem.hh>
//#include<dune/pdelab/newton/newton.hh>

#include <dune/grid/io/file/gmshreader.hh>
#include <dune/grid/io/file/vtk/subsamplingvtkwriter.hh>

#include "bctype.hh"
#include "operator_simple_iter.hh"
#include "l2norm.hh"

	
template<typename GV>
void driver_simple_iter(const GV& gv, double H){
	using namespace Dune::PDELab;  // Da skratimo imena

	// Prostor konačnih elemenata (grid function space) i GridOperator.
	const int dim = GV::dimension;
	const int k = 1;
	typedef QkLocalFiniteElementMap<GV,double,double,k>               FEM;
	typedef ConformingDirichletConstraints                            CONSTRAINTS;
	typedef ISTL::VectorBackend<>                                     VBE;
	typedef GridFunctionSpace<GV,FEM,CONSTRAINTS,VBE>                 GFS;
	typedef typename GFS::template ConstraintsContainer<double>::Type CC;
	typedef Backend::Vector<GFS,double>                               U;  // Na ovaj način U definiramo prije GO
	typedef DiscreteGridFunction<GFS,U>                               DGF;
	typedef DiscreteGridFunctionGradient<GFS,U>                       DGG;
	typedef SeqMinimalSurfaceLOP<DGG, FEM>                            LOP;
	typedef ISTL::BCRSMatrixBackend<>                                 MBE;
	typedef GridOperator<
		GFS,GFS,              /* prostor KE rješenja i test funkcije */
		LOP,                  /* lokalni operator */
		MBE,                  /* matrix backend */
		double,double,double, /* tipovi u domeni, slici i jakobijanu */
		CC,CC                 /* ograničenja za prostor rješenja i test funkcija. */
		> GO;

	FEM fem(gv);
	GFS gfs(gv,fem);
	CC cc;
	DirichletBdry bctype;
	constraints(bctype, gfs, cc); // asembliranje ograničenja Dirichletovog tipa
	std::cout << "constrained dofs=" << cc.size() << " of " << gfs.globalSize() << std::endl;
	U u(gfs,0.0), u_prev(gfs,0.0);
	DGF udgf(gfs,u), u_prevdgf(gfs,u_prev);
	DGG gradu_prevdgf(gfs,u_prev);
	LOP lop(H, gradu_prevdgf, fem);  // LOP uzima gradijent prethodnog rješenja
	MBE mbe(9);  // traži prosječan broj ne-nul elemenata u redu (=9)
	GO go(gfs,cc,gfs,cc,lop,mbe);

	//  Konstrukcija rješavača.
	typedef BCExtension<GV,double>                 G;
	typedef ISTLBackend_SEQ_BCGS_SSOR              LS;
	typedef StationaryLinearProblemSolver<GO,LS,U> SLP;

	G g(gv);
	interpolate(g, gfs, u);
	interpolate(g, gfs, u_prev);
	LS ls(5000,true);        // max 5000 iteracija, verbosity = true
	SLP slp(go,ls,u,1e-10);  // redukcija = 1e-10

	int i = 0;
	while(i < 500){
		std::cout << std::endl;
		std::cout << "Step " << i << std::endl;
		slp.apply();
		if (l2normdiff(gv, udgf, u_prevdgf) / l2norm(gv, u_prevdgf) < 1e-10) break;
		u_prev = u;
		i++;
	}
	double error = l2normdiff(gv, udgf, u_prevdgf) / l2norm(gv, u_prevdgf);
	std::cout << "Iterations: " << i 
			  << "\nError: " << error << std::endl;

	// Grafički izlaz (VTK)
	Dune::VTKWriter<GV> vtkwriter(gv,Dune::VTK::conforming);
	vtkwriter.addVertexData(std::make_shared<VTKGridFunctionAdapter<DGF>>(udgf,"seq_elevation"));
	vtkwriter.write("seq_surface", Dune::VTK::ascii); //Dune::VTK::appendedraw);

}

#endif

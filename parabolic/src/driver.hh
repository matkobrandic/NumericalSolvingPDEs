#ifndef DRIVER_HH_17294361
#define DRIVER_HH_17294361
#include <dune/pdelab/finiteelementmap/qkfem.hh>
#include <dune/pdelab/backend/istl.hh>
#include <dune/pdelab/stationary/linearproblem.hh>
#include <dune/pdelab/instationary/onestep.hh>
#include <dune/istl/bvector.hh>
#include <dune/istl/operators.hh>
#include <dune/istl/solvers.hh>
#include <dune/istl/preconditioners.hh>
#include <dune/istl/io.hh>
#include <dune/istl/superlu.hh>
#include <dune/pdelab/common/function.hh>
#include <dune/pdelab/common/vtkexport.hh>
#include <dune/pdelab/constraints/conforming.hh>
#include <dune/pdelab/gridfunctionspace/gridfunctionspace.hh>
#include <dune/pdelab/gridfunctionspace/gridfunctionspaceutilities.hh>
#include <dune/pdelab/gridfunctionspace/interpolate.hh>
#include<dune/pdelab/newton/newton.hh>
#include <dune/grid/io/file/vtk/vtksequencewriter.hh>
#include <dune/pdelab/gridoperator/gridoperator.hh>
#include <dune/pdelab/gridoperator/onestep.hh>
#include "bctype.hh"
#include "bcextension.hh"
#include "space_operator.hh"
#include "time_operator.hh"

/** Upravljačka rutina koja koordinira sav posao osim konstrukciju 
 *  mreže. 
 *  @tparam GV = Leaf grid view tip 
 *
 *  @param gv = leaf grid view 
 *  @param dt = početni vremenski korak
 *  @param tend = vrijeme simulacije
 *  */
template<class GV>
void driver(const GV& gv, double dt, double tend){
	const int Qk_order = 2;
	using FEM = Dune::PDELab::QkLocalFiniteElementMap<GV,double,double,Qk_order>;
	using CON = Dune::PDELab::ConformingDirichletConstraints;
	using VBE = Dune::PDELab::ISTL::VectorBackend<>;
	using GFS = Dune::PDELab::GridFunctionSpace<GV,FEM,CON,VBE>;
	using  CC = typename GFS::template ConstraintsContainer<double>::Type;
	// Vektor stupnjeva slobode
	using U =  Dune::PDELab::Backend::Vector<GFS, double>;
	// Lokalni operator za prostorni dio
	using SLOP = StationaryLocalOperator<BCTypeParam,FEM>;
	// Lokalni operator za vremenski dio
	using TLOP = TimeLocalOperator<FEM>;
	using MBE = Dune::PDELab::ISTL::BCRSMatrixBackend<>;
	// Grid operatori koji odgovaraju vremenskom i prostornom dijelu
	using GO0 = Dune::PDELab::GridOperator<GFS,GFS,SLOP,MBE,double,double,double,CC,CC>;
	using GO1 = Dune::PDELab::GridOperator<GFS,GFS,TLOP,MBE,double,double,double,CC,CC>;
	// Potpuni grid operator se dobiva jednokoračnom vremenskom metodom
	using IGO = Dune::PDELab::OneStepGridOperator<GO0,GO1> ;
	// Vektor stupnjeva slobode može se dobiti i na ovaj način:
	//  using U = typename IGO::Traits::Domain;
	using G = BCExtension<GV>;      // rubni uvjet i početni uvjet
	// Linearni solver
	using LS = Dune::PDELab::ISTLBackend_SEQ_BCGS_SSOR;
	// Zadaća je nelinearna
	using PDESOLVER = Dune::PDELab::Newton<IGO,LS,U>;
	// Tip vremenske diskretizacije
	using OSM = Dune::PDELab::OneStepMethod<double,IGO,PDESOLVER,U,U>;
	// Parametri vremenske diskretizacije
	using TDM = Dune::PDELab::Alexander2Parameter<double>;
	// Druge mogućnosti:
	// using TDM = Dune::PDELab::ImplicitEulerParameter<double>; // implicitni Euler
	//  using TDM = Dune::PDELab::OneStepThetaParameter<double>; // theta metoda, theta parametar uzima konstruktor
	using DGF = Dune::PDELab::DiscreteGridFunction<GFS,U>;
	// Za Q2 elemente trebamo koristiti subsampling.
	using VTKW = Dune::SubsamplingVTKWriter<GV>;
	//  using VTKW = Dune::VTKWriter<GV>;
	double time = 1.0; // početni trenutak

	FEM fem(gv);
	GFS gfs(gv,fem);
	BCTypeParam bctype;     // Korisnička klasa koja detektora Dirichletov rubni uvjet.
	bctype.setTime(time);   // b.c. ovisi o vremenu
	CC cc;
	Dune::PDELab::constraints( bctype, gfs, cc ); // Definiraj Dirichletove vrhove

	SLOP slop(bctype,4);        // prostorni lokalni operator
	TLOP tlop(4);               // vremenski lokalni operator

	MBE mbe(9);
	GO0 go0(gfs, cc, gfs, cc, slop, mbe);  // prostorni GO
	GO1 go1(gfs, cc, gfs, cc, tlop, mbe);  // vremenski GO
	IGO igo(go0,go1);           // grid operator jednokoračne metode

	// typename GO0::Traits::Jacobian jac(go0);
	// std::cout << jac.patternStatistics() << std::endl;

	U uold(gfs,0.0);             // rješenje u t=t^n, početni uvjet je 0
	// U exact(gfs,0.0); 
	G bcond(gv);                 // rubni uvjet sada ovisi o vremenu
	bcond.setTime(time);
	Dune::PDELab::interpolate(bcond,gfs,uold);  // početni uvjet je dan ovdje
	//Dune::PDELab::copy_constrained_dofs(cc,exact,uold);

	// egzaktno rješenje
	U exact(gfs, 0.0);
	ExactGF<GV> exact_dgf(gv);
	exact_dgf.setTime(time);
	Dune::PDELab::interpolate(exact_dgf, gfs, exact);

	LS ls(5000, false);
	PDESOLVER pdesolver(igo,uold,ls);
	pdesolver.setReassembleThreshold(0.0);
	pdesolver.setVerbosityLevel(1);
	pdesolver.setReduction(1e-10);
	pdesolver.setMinLinearReduction(1e-4);
	pdesolver.setMaxIterations(25);
	pdesolver.setLineSearchMaxIterations(10);

	TDM method;                      // koeficijenti vremenske diskretizacije
	OSM osm(method,igo,pdesolver);   // jednokoračna metoda za rješenje sustava
	osm.setVerbosityLevel(1);

	// Iscrtavanje
	DGF udgf(gfs,uold);
	VTKW vtkwriter(gv, Dune::RefinementIntervals{2});
	vtkwriter.addVertexData(std::make_shared<Dune::PDELab::VTKGridFunctionAdapter<DGF>>(udgf,"solution"));
	vtkwriter.addVertexData(std::make_shared<Dune::PDELab::VTKGridFunctionAdapter<ExactGF<GV>>>(exact_dgf, "exact"));
	Dune::VTKSequenceWriter<GV> writer(std::make_shared<VTKW>(vtkwriter), "out");
	writer.write(1.0);

	// Vremenska petlja
	U unew(gfs,0.0);        // sljedeći vremenski sloj -jednokoračna metoda

	while(time < tend){
		// postavi novo vrijeme u BC klasu -- za Dirichletovu granicu koja se mijenja u vremenu
		//  bctype.setTime(time+dt);                   // izračunaj Dirichletov r.u
		//  cc.clear();                                // u ovom vremenskom trenutku
		//  Dune::PDELab::constraints(bctype,gfs,cc);

		// Riješi sustav. Metoda apply će interpolirati i rubni uvjet
		osm.apply(time, dt, uold, bcond, unew);  // novo rješenje je u vektru unew
		// Broj linearnih iteracija
		int noIter = osm.getPDESolver().result().iterations;   
		std::cout << "  === No of linear iterations = " << noIter << ", dt = " << dt << std::endl;
		// Pripremi sljedeći vremenski korak
		uold = unew;
		time += dt;
		exact_dgf.setTime(time);    
		writer.write(time);
	}
}

#endif

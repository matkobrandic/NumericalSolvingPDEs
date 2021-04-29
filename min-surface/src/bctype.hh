#ifndef BCTYPE_HH
#define BCTYPE_HH

#include <dune/common/fvector.hh>
#include <dune/pdelab/common/function.hh>
#include <dune/pdelab/constraints/common/constraintsparameters.hh>


// Rubni uvjet. Ovdje je zadan kao primjer za testiranje. Kada je zadatak gotov
// zamijeniti ga nekim svojim primjerom.
template <int dim>
double bc(Dune::FieldVector<double, dim> const & glob){
	double x = glob[0], y=glob[1];
	return std::sin(x)*std::sin(x) + std::cos(y)*std::cos(y);
}



// Klasa za identifikaciju Dirichetove granice. Ne treba mijenjati.
class DirichletBdry : public Dune::PDELab::DirichletConstraintsParameters{
public:
	//  intersection = stranica elementa (u 3D) ili brid elementa (u 2D)
	//  coord        = lokalne koordinate točke na "intersectionu" koja se ispituje
	//  povratna vrijednost: true ako je točka na Dirichletovoj granici
	//                       false ako nije. 
	template<typename I>
	bool isDirichlet(const I & intersection,
					 const Dune::FieldVector<typename I::ctype, I::mydimension> & coord  // 2.6.0 mydimension daje dim interfacea
					) const{
	//Dune::FieldVector<typename I::ctype, I::coorddimension> xg = intersection.geometry().global( coord );
	return true;  // Dirichletov uvjet na cijeloj granici
	}
};

// Klasa koja daje Dirichletov rubni uvjet. Ne treba mijenjati.
template<typename GV, typename RF>
class BCExtension
	: public Dune::PDELab::GridFunctionBase<Dune::PDELab::
										  GridFunctionTraits<GV,RF,1,Dune::FieldVector<RF,1> >, 
										  BCExtension<GV,RF> > {
	const GV& gv;
public:
	typedef Dune::PDELab::GridFunctionTraits<GV,RF,1,Dune::FieldVector<RF,1> > Traits;

	// Konstruktor samo uzima referencu na  GridView objekt. 
	BCExtension (const GV& gv_) : gv(gv_) {}

	// Izračunaj Dirichletovu vrijednost na elementu. Ako točka nije na 
	// Dirichletovoj granici, nda funkcija daje proširenje Dirichletovog rubnog
	// uvjeta na čitavu domenu. To je proširenje u osnovi proizvoljno. 
	// e     = element 
	// xlocal = lokalne koordinate točke u kojoj se računa Dirichletova vrijednost
	// y      = izračunata Dirichletova vrijednost
	inline void evaluate (const typename Traits::ElementType& e,
						  const typename Traits::DomainType& xlocal,
						  typename Traits::RangeType& y) const{
		const int dim = Traits::GridViewType::Grid::dimension;
		typedef typename Traits::GridViewType::Grid::ctype ctype;

		// Pretvori lokalne koordinate u globalne
		Dune::FieldVector<ctype,dim> x = e.geometry().global(xlocal);
		y = bc(x);
	}

	// Vrati referencu na GridView
	inline const GV& getGridView () {return gv;}
};

#endif

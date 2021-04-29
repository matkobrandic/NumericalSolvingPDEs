#include <dune/pdelab/constraints/conforming.hh>

/*  Selekcija Dirichletove granice */
class BCTypeParam
  : public Dune::PDELab::DirichletConstraintsParameters {
	// Dirichletova granica se može micati u vremenu.
	double time;
public:
	template<typename I>
	bool isDirichlet(const I & intersection,   
	const Dune::FieldVector<typename I::ctype, I::coorddimension-1> & coord) const{
    	return true; // sve ostalo je
	}
	//! Postavi vrijeme.
	void setTime (double t) { time = t; }
};

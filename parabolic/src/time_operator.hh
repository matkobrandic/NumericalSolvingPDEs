#include <dune/geometry/referenceelements.hh>
#include <dune/geometry/quadraturerules.hh>
#include <dune/pdelab/common/geometrywrapper.hh>
#include <dune/pdelab/localoperator/defaultimp.hh>
#include <dune/pdelab/localoperator/pattern.hh>
#include <dune/pdelab/localoperator/flags.hh>
#include <dune/pdelab/localoperator/idefault.hh>
#include <dune/pdelab/finiteelement/localbasiscache.hh>

/** Vremenski lokalni operator. Ovdje je to samo skalarni produkt.
 *
 * \f{align*}{
                \int_\Omega uv dx
 * \f}
 */
template <typename FEM>
class TimeLocalOperator 
	:	public Dune::PDELab::NumericalJacobianApplyVolume<TimeLocalOperator<FEM>>,
		public Dune::PDELab::NumericalJacobianVolume<TimeLocalOperator<FEM>>,
		public Dune::PDELab::FullVolumePattern,
		public Dune::PDELab::LocalOperatorDefaultFlags,
		public Dune::PDELab::InstationaryLocalOperatorDefaultMethods<double>{
public:
	enum { doPatternVolume = true };
	enum { doAlphaVolume = true };

	TimeLocalOperator (unsigned int intorder_=2)
		: intorder(intorder_), time(0.0)
	{}

	//! set time for subsequent evaluation
	void setTime(double t){
		time = t;
	}

	// volume integral depending on test and ansatz functions
	template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
	void alpha_volume(const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, R& r) const{
		const int dim = EG::Geometry::mydimension;
		// select quadrature rule
		Dune::GeometryType gt = eg.geometry().type();
		const auto& rule = Dune::QuadratureRules<double,dim>::rule(gt,intorder);
		// loop over quadrature points
		for(auto const & ip : rule){
			// Bazne funkcije u integracijskoj točki.
			auto& phi = cache.evaluateFunction(ip.position(),lfsu.finiteElement().localBasis());
			// Rješenje u
			double u = 0.0;
			for(size_t i = 0; i < lfsu.size(); ++i){
				u += x(lfsu,i)*phi[i];
			}
		    // u*phi_i
		    double factor = ip.weight() * eg.geometry().integrationElement(ip.position());
		    for(size_t i = 0; i < lfsv.size(); ++i){
				r.accumulate(lfsv,i,u*phi[i]*factor);
		    }
		}
	}
private:
	unsigned int intorder;
	double time;
	typedef typename FEM::Traits::FiniteElementType::Traits::LocalBasisType LocalBasis;
	Dune::PDELab::LocalBasisCache<LocalBasis> cache;
};

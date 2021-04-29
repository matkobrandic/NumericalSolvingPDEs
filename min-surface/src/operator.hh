#ifndef OPERATOR_HH
#define OPERATOR_HH

#include <dune/geometry/quadraturerules.hh>
#include <dune/geometry/referenceelements.hh>

#include <dune/pdelab/localoperator/defaultimp.hh>
#include <dune/pdelab/localoperator/flags.hh>
#include <dune/pdelab/localoperator/pattern.hh>
#include <dune/pdelab/finiteelement/localbasiscache.hh>

//#include "bctype.hh"
/** Lokalni operator za zadaću :
 *
 *    div(  grad u / (1+ |grad u|^2)^(1/2) )  = 2H   u \Omega
 *                   u = g(x)   na \partial\Omega
 *
 * sa konformnim konačnim elementima svih tipova u  2D.
 *
 * \tparam BCType klasa koja indicira rubni uvjet
 */



template<typename BCType, typename FEM>
class MinimalSurfaceLOP : // derivacijska lista -- jakobijan i pattern računa PDELab
  public Dune::PDELab::NumericalJacobianApplyVolume  <MinimalSurfaceLOP<BCType, FEM> >,
  public Dune::PDELab::NumericalJacobianVolume       <MinimalSurfaceLOP<BCType, FEM> >,
  public Dune::PDELab::NumericalJacobianApplyBoundary<MinimalSurfaceLOP<BCType, FEM> >,
  public Dune::PDELab::NumericalJacobianBoundary     <MinimalSurfaceLOP<BCType, FEM> >,
  public Dune::PDELab::FullVolumePattern,
  public Dune::PDELab::LocalOperatorDefaultFlags
{
public:
	// Zastavice koje signaliziraju da na svakom elementu treba zvati:
	enum { doPatternVolume = true };  // metodu za računanje patterna (iz volumnih doprinosa)
	enum { doAlphaVolume = true };    // alpha_volume
	using  LocalBasis = typename FEM::Traits::FiniteElementType::Traits::LocalBasisType;

	MinimalSurfaceLOP(double H_, const BCType& bctype_, // boundary cond.type
					  const FEM & fem_,
					  unsigned int intorder_=3) :
		H(H_), bctype( bctype_ ), fem(fem_), intorder( intorder_ )
	{}

	// Računanje volumnog integrala
	// eg   = element (geometry)
	// lfsu = lokalni prostor funkcija za rješenje
	// lfsv = lokalni prostor funkcija za test funkciju
	// x    = vektor koeficijenata rješenja
	// r    = lokalni rezidual
	template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
	void alpha_volume (const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, R& r) const{
		const int dim= EG::Geometry::mydimension;
		using Gradient = Dune::FieldVector<double, dim>;

		auto gt = eg.geometry().type();
		auto & rule = Dune::QuadratureRules<double, dim>::rule(gt, intorder);

		for(auto qpoint : rule){
			auto const & xi = qpoint.position();
			auto & phihat = cache.evaluateFunction(xi, lfsu.finiteElement().localBasis());
			auto & gradphihat = cache.evaluateJacobian(xi, lfsu.finiteElement().localBasis());

			double u = 0.0;
			for(int i=0; i<lfsu.size(); ++i){
				u += x(lfsu, i)*phihat[i];
			}

			auto const & jac = eg.geometry().jacobianInverseTransposed(xi);
			std::vector<Gradient> gradphi(lfsu.size());
			for(int i=0; i<lfsu.size(); ++i){
				jac.mv(gradphihat[i][0], gradphi[i]);
			}

			Gradient grad_u(0.0);
			for(int i=0; i<lfsu.size(); ++i){
				grad_u.axpy(x(lfsu, i), gradphi[i]); // grad_u += x(lfsu, i) * gradphi[i]
			}

			auto dx = qpoint.weight() * eg.geometry().integrationElement(xi);
			for(int i=0; i<lfsu.size(); ++i){
				r.accumulate(lfsu, i, ( - (grad_u*gradphi[i]) / (std::sqrt(1 + grad_u.dot(grad_u))) - 2*H*phihat[i] )*dx );
			}
		}
	}

private:
	const double H; // srednja zakrivlenost plohe
	BCType const & bctype;
	FEM const & fem;
	unsigned int intorder;
	Dune::PDELab::LocalBasisCache<LocalBasis> cache;
};
#endif

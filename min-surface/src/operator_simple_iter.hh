#ifndef OPERATOR_SIMPLE_ITER_HH
#define OPERATOR_SIMPLE_ITER_HH

#include <dune/geometry/quadraturerules.hh>
#include <dune/geometry/referenceelements.hh>

#include <dune/pdelab/localoperator/defaultimp.hh>
#include <dune/pdelab/localoperator/flags.hh>
#include <dune/pdelab/localoperator/pattern.hh>
#include <dune/pdelab/finiteelement/localbasiscache.hh>

//#include "bctype.hh"
/** Lokalni operator za zadaću :
 *
 *   - div(  grad u / (1+ |grad u|^2)^(1/2) )  = 0   u \Omega
 *                   u = g(x)   na \partial\Omega
 *
 * sa konformnim konačnim elementima svih tipova u  2D.
 *
 * \tparam DGG klasa koja indicira rubni uvjet
 */



template<typename DGG, typename FEM>
class SeqMinimalSurfaceLOP : // derivacijska lista -- jakobijan i pattern računa PDELab
	public Dune::PDELab::NumericalJacobianApplyVolume  <SeqMinimalSurfaceLOP<DGG, FEM> >,
	public Dune::PDELab::NumericalJacobianVolume       <SeqMinimalSurfaceLOP<DGG, FEM> >,
	public Dune::PDELab::NumericalJacobianApplyBoundary<SeqMinimalSurfaceLOP<DGG, FEM> >,
	public Dune::PDELab::NumericalJacobianBoundary     <SeqMinimalSurfaceLOP<DGG, FEM> >,
	public Dune::PDELab::FullVolumePattern,
	public Dune::PDELab::LocalOperatorDefaultFlags
	{
	public:
		// Zastavice koje signaliziraju da na svakom elementu treba zvati:
		enum { doPatternVolume = true };  // metodu za računanje patterna (iz volumnih doprinosa)
		enum { doAlphaVolume = true };    // alpha_volume
		using  LocalBasis = typename FEM::Traits::FiniteElementType::Traits::LocalBasisType;

		SeqMinimalSurfaceLOP(double H_, const DGG& dgfgrad_, // boundary cond.type
							const FEM & fem_,
							unsigned int intorder_=3) :
			H(H_), dgfgrad( dgfgrad_ ), fem(fem_), intorder( intorder_ )
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
			typename Dune::PDELab::FunctionTraits<double, dim, Dune::FieldVector<double, dim>, double, dim, Dune::FieldVector<double, dim> >::RangeType grad_u_prev;


			auto gt = eg.geometry().type();
			auto & rule = Dune::QuadratureRules<double, dim>::rule(gt, intorder);

			for(auto qpoint : rule){
				auto const & xi = qpoint.position();
				auto & phihat = cache.evaluateFunction(xi, lfsu.finiteElement().localBasis());
				auto & gradphihat = cache.evaluateJacobian(xi, lfsu.finiteElement().localBasis());

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
				dgfgrad.evaluate(eg.entity(), xi, grad_u_prev); // dgfgrad primam, računam njegovu vrijednost u točki integracije

				for(int i=0; i<lfsu.size(); ++i){
					r.accumulate(lfsu, i, ( - (grad_u*gradphi[i]) / (std::sqrt(1 + grad_u_prev.dot(grad_u_prev))) - 2*H*phihat[i] )*dx );
				}
			}
		}

	private:
		const double H;
		DGG const & dgfgrad;
		FEM const & fem;
		unsigned int intorder;
		Dune::PDELab::LocalBasisCache<LocalBasis> cache;
};
#endif

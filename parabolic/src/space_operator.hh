#include<dune/pdelab/localoperator/idefault.hh>
#include<dune/geometry/quadraturerules.hh>
#include<dune/geometry/referenceelements.hh>
#include<dune/pdelab/localoperator/defaultimp.hh>
#include<dune/pdelab/localoperator/flags.hh>
#include<dune/pdelab/localoperator/pattern.hh>
#include<dune/pdelab/finiteelement/localbasiscache.hh>


/** Lokalni operator za zadaću:
 *
 *   - \Delta u + a*u = f   u \Omega
 *                  u = g   na \Gamma_D\subseteq\partial\Omega
 *  -\nabla u \cdot n = j   na \Gamma_N = \partial\Omega\setminus\Gamma_D
 *
 *
 * \tparam BCType - klasa koja određuje tip rubnog uvjeta.
 */
template <typename BCType, typename FEM>
class StationaryLocalOperator :
	public Dune::PDELab::NumericalJacobianApplyVolume<StationaryLocalOperator<BCType,FEM> >,
	public Dune::PDELab::NumericalJacobianVolume<StationaryLocalOperator<BCType,FEM> >,
	public Dune::PDELab::NumericalJacobianApplyBoundary<StationaryLocalOperator<BCType,FEM> >,
	public Dune::PDELab::NumericalJacobianBoundary<StationaryLocalOperator<BCType,FEM> >,
	public Dune::PDELab::FullVolumePattern,
	public Dune::PDELab::LocalOperatorDefaultFlags,
	public Dune::PDELab::InstationaryLocalOperatorDefaultMethods<double>
{
public:
	enum { doPatternVolume = true };
	enum { doAlphaVolume = true };
	enum { doAlphaBoundary = false };

	StationaryLocalOperator(BCType& bctype_, // boundary cond.type
							unsigned int intorder_=2) :
							bctype( bctype_ ), intorder( intorder_ )
  {}

	// volume integral depending on test and ansatz functions
	template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
	void alpha_volume(const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, R& r) const{
		//lfsu == lfsv
		const int dim = EG::Geometry::coorddimension;
		using Gradient = Dune::FieldVector<double,dim>;

		auto gt = eg.geometry().type();
		const auto & rule = Dune::QuadratureRules<double,dim>::rule(gt, intorder);

		for(const auto & ip : rule){
			// Izračunaj bazne funkcije
			auto& phi = cache.evaluateFunction(ip.position(), lfsu.finiteElement().localBasis());
			// rješenje
			double u = 0.0;
			for (size_t i = 0; i < lfsu.size(); ++i){
				u += x(lfsu,i) * phi[i];
			}
			// Gradijenti baznih funkcija
			auto& gradphihat = cache.evaluateJacobian(ip.position(), lfsu.finiteElement().localBasis());
			const auto & jac = eg.geometry().jacobianInverseTransposed(ip.position());
			// Gradijenti baznih funkcija na fizičkom elementu.
			std::vector<Gradient> gradphi(lfsu.size());
			for(size_t i=0; i<lfsu.size(); i++){
				jac.mv(gradphihat[i][0],gradphi[i]);
			}
			// Gradijent rješenja u
			Gradient gradu(0.0);
			for(size_t i=0; i<lfsu.size(); ++i){
				gradu.axpy(x(lfsu,i),gradphi[i]);
			}
			double m = 2;
			// integrate grad u * grad phi_i + a*u*phi_i - f phi_i
			double factor = ip.weight()*eg.geometry().integrationElement(ip.position());
			for (size_t i = 0; i < lfsv.size(); i++)
			r.accumulate(lfsv, i, (pow(u, m - 1) * m * (gradu * gradphi[i]) ) * factor);
		}
	}

	void preStep(double time, double dt, int stages){
    	bctype.setTime(time+dt); // ako treba promijeni Dirichletovu granicu
    	Dune::PDELab::InstationaryLocalOperatorDefaultMethods<double>::preStep(time, dt, stages);
	}
private:
	BCType& bctype;
	unsigned int intorder;
	typedef typename FEM::Traits::FiniteElementType::Traits::LocalBasisType LocalBasis;
	Dune::PDELab::LocalBasisCache<LocalBasis> cache;
};

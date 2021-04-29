#include <cmath>
#include <dune/pdelab/common/function.hh>

template <int dim>
double exact(Dune::FieldVector<double, dim> const & glob, double time){
	double m = 2;
	double alfa = dim / (dim * (m-1) + 2);
	double k = (alfa * (m-1)) / (2 * m * dim);
	double beta = alfa / dim;
	double temp1 = 1 / (pow(time, alfa));
	double norma = glob.two_norm2();
	double C = 10000 * k / pow(time, 2 * beta);
	double temp2 = C - k * (norma) / (pow(time, 2 * beta));
	temp2 = std::max(temp2, 0.0);
	temp2 = temp1 * pow(temp2, 1/(m - 1));
	return temp2;
}

// Klasa obilježja. Ovo je samo pokrata koja čini kod čitljivijim.
template<typename GV>
using ScalarTraits = Dune::PDELab::GridFunctionTraits<GV,double,1,Dune::FieldVector<double,1>>;


/*  Dirichletov rubni uvjet proširen na čitavu domenu.
 *  U t=0 daje inicijalni uvjet!
 */
template<typename GV>
class BCExtension
  : public Dune::PDELab::GridFunctionBase<ScalarTraits<GV>, BCExtension<GV> >{
	const GV& gv;
	double  time;
	public:
	using Traits = ScalarTraits<GV>;

	// Sačuvaj gridview
	BCExtension (const GV& gv_) : gv(gv_) {}

	inline void evaluate (const typename Traits::ElementType& e, 
						  const typename Traits::DomainType& xlocal,
						  typename Traits::RangeType& y) const{
		auto x = e.geometry().global(xlocal);
		if (time == 1){
			y = exact(x,time);
		}
		else{
			y = 0.0;
		}
		return;
	}

// Referenca na gridview
 inline const GV& getGridView () {return gv;}

void setTime (double t) {
	time = t;}
 }; 
template <typename GV>
using ATraits = Dune::PDELab::AnalyticGridFunctionTraits<GV,double,1>;// 1 = skalarna funkcija
template <typename GV>

class ExactGF : public Dune::PDELab::AnalyticGridFunctionBase<ATraits<GV>, ExactGF<GV>>{
	double time;
	public:typedef Dune::PDELab::AnalyticGridFunctionBase<ATraits<GV>, ExactGF<GV> > BaseT;
	ExactGF(GV const & gv) : BaseT(gv) {}

	void evaluateGlobal(typename ATraits<GV>::DomainType const & x, typename ATraits<GV>::RangeType & y) const{
		y = exact(x, time);
	}
	// Postavljanje vremena. Potrebno pozvati prije evaluate() !
	void setTime (double t){
		time = t;
	}
};

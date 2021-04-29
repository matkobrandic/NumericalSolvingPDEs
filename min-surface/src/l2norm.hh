#ifndef L2NORM_HH_IS_INCLUDED
#define L2NORM_HH_IS_INCLUDED
#include <dune/geometry/quadraturerules.hh>

// L2 norma razlige dgf1 i dgf2
template <typename GV, typename DGF>
double l2normdiff(GV const & gv, DGF const & dgf1,  DGF const & dgf2, int p = 3){
	double value = 0.0;
	const int dim = GV::dimension;
	for(auto const & element : elements(gv)){
		typename Dune::PDELab::FunctionTraits<double, dim, Dune::FieldVector<double, dim>, double, 1, Dune::FieldVector<double, 1> >::RangeType v1;
		typename Dune::PDELab::FunctionTraits<double, dim, Dune::FieldVector<double, dim>, double, 1, Dune::FieldVector<double, 1> >::RangeType v2;
		const auto geometry = element.geometry();
		const auto gt = geometry.type();
		const auto& rule = Dune::QuadratureRules<double,dim>::rule(gt,p);
		if (rule.order() < p) std::exit(-1);

		double res = 0.0;
		for(auto const & qpoint : rule){
			dgf1.evaluate(element, qpoint.position(), v1);
			dgf2.evaluate(element, qpoint.position(), v2);
			double weight = qpoint.weight();
			double detjac = geometry.integrationElement(qpoint.position());
			res += (v1-v2) * (v1-v2) * weight * detjac;
		}
		value += res;
	}
	return std::sqrt(value);
}

// L2 norma
template <typename GV, typename DGF>
double l2norm(GV const & gv, DGF const & dgf1, int p = 3){
	double value = 0.0;
	const int dim = GV::dimension;
	for(auto const & element : elements(gv)){
		typename Dune::PDELab::FunctionTraits<double, dim, Dune::FieldVector<double, dim>, double, 1, Dune::FieldVector<double, 1> >::RangeType v1;
		const auto geometry = element.geometry();
		const auto gt = geometry.type();
		const auto& rule = Dune::QuadratureRules<double,dim>::rule(gt,p);
		if (rule.order() < p) std::exit(-1);

		double rezultat = 0.0;
		for(auto const & qpoint : rule){
			dgf1.evaluate(element, qpoint.position(), v1);
			double weight = qpoint.weight();
			double detjac = geometry.integrationElement(qpoint.position());
			rezultat += v1 * v1 * weight * detjac;
		}
		value += rezultat;
	}
	return std::sqrt(value);
}

#endif


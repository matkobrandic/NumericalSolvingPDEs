#ifndef DUNE_GRID_HOWTO_TRANSPORTPROBLEM2_HH
#define DUNE_GRID_HOWTO_TRANSPORTPROBLEM2_HH

#include <dune/common/fvector.hh>
//===============================================================
// Zadaća:
//
//   dc/dt + div( c q(x,t) ) = 0   u Omega
//                   c = b   na Gamma_in
//                   c = c0  za t=0
//
// Definiramo funkcije c0(x), q(x,t) i b(x,t).
//===============================================================

// OVDJE NE TREBA NIŠTA MIJENJATI

// početna koncentracija
template<int dimworld>
double c0 (const Dune::FieldVector<double, dimworld>& x)
{
  if (x.two_norm()>0.125 && x.two_norm()<0.5)
    return 1.0;
  else
    return 0.0;
}

// Rubni uvjet na ulaznoj granici.
template<int dimworld>
double b(const Dune::FieldVector<double,dimworld>& x, double t)
{
  return 0.0;
}

// vektor brzine q(x,t)
template<int dimworld>
Dune::FieldVector<double,dimworld> q(const Dune::FieldVector<double,dimworld>& x, double t)
{
  Dune::FieldVector<double,dimworld> q(0.0);
  q[0] = 1.0;
  return q;
}
#endif // DUNE_GRID_HOWTO_TRANSPORTPROBLEM2_HH

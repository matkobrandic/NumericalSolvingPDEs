#ifndef DUNE_GRID_HOWTO_EVOLVE_HH
#define DUNE_GRID_HOWTO_EVOLVE_HH

#include <dune/common/fvector.hh>

/*  Rutina koja računa evoluciju koncentracije "c" s vremenskog sloja t na vremenski sloj 
 *  t+dt. Korak dt se izračunava na osnovu CFL uvjeta i vraća u glavni program.
 *
 *  gv = leaf grid view
 *  idxSet = odgovarajući indexSet
 *  c = na ulazu: vektor rješenja na prethodnom vremenskom  sloju
 *      na izlazu: vektor rješenja na sljedećem vremenskom sloju.
 *  t = trenutno vrijeme
 *  dt = vremenski korak je izlazna vrijednost. Ulazna vrijednost se ignorira.
 *
 */
template<class GV, class IS, class V>
void evolve (GV const & gv, const IS& idxSet, V& c, double t, double& dt)
{
  // Treba nam privremeni vektor za updatiranje jer dt računamo istovremeno.
  // novi_c = stari_c + dt*update
  V update(c.size());                                  
  for (typename V::size_type i=0; i<c.size(); i++) update[i] = 0;


  //  Petlja po svim elementima
  for(auto const & element : elements(gv))
  {
      // VAŠ KOD DOLAZI OVDJE. 
      // Ovdje treba izračunati vektor update[] i dt određen sa CFL uvjetom
      throw std::runtime_error("Kod još nije napisan."); 
  } 

  // dodajemo izračunatu korekciju na vektor c.
  for (unsigned int i=0; i<c.size(); ++i)
    c[i] += dt*update[i];                         

  return;
}

#endif //DUNE_GRID_HOWTO_EVOLVE_HH

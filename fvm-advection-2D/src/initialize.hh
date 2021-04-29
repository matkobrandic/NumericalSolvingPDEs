#ifndef DUNE_GRID_HOWTO_INITIALIZE_HH
#define DUNE_GRID_HOWTO_INITIALIZE_HH

#include <dune/common/fvector.hh>

/**
 *    Početni uvjet.
 *    Inicijaliziraj vektor c s početnim uvjetom
 *    GV = GridView
 *    IS = IndexSet
 *    V = tip vektora
 *
 */
template<class GV, class IS, class V>
void initialize (GV& gv, const IS& idxSet, V& c)
{
   // VAŠ KOD DOLAZI OVDJE
   // Ovdje staviti kod za inicijalizaciju vektora c.
   throw std::runtime_error("Kod još nije napisan");  
}

#endif // DUNE_GRID_HOWTO_INITIALIZE_HH

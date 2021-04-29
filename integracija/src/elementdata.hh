#ifndef ELEMENTDATA_HH_
#define ELEMENTDATA_HH_

#include <dune/grid/common/mcmgmapper.hh>
#include <vector>
// Ova funkcija uzima grid i funkciju zadanu kao funkcijski objekt te formira vektor
// coeff koji se sastoji od vrijednosti funkcije u centrima elemenata.
// Na primjer, ako je x_i centar elementa K_i, onda je
//                   coeff[i] = f(x_i).
//
// LGV = LeafGridView
// F = Funkcijski objekt
template<class LGV, class F>
void elementdata (const LGV& gridView, const F& f, std::vector<double> & coeff)
{
  // Kreiramo mapper koji indeksira sve elemente mreže (bez obzira na njihov tip).
  // mapper treba referencu na Grid (a ne GridView).
  //Dune::LeafMultipleCodimMultipleGeomTypeMapper<typename LGV::Grid, Dune::MCMGElementLayout> mapper(gridView.grid());
    //KORISTIMO set UMJESTO mapper pa ih svugdje zamijenimo
  auto const & set = gridView.indexSet();
  // Dajmo vektoru koeficijenata odgovarajuću dimenziju
  //coeff.resize(mapper.size());
  coeff.resize(set.size(0));

  // Iteriramo kroz sve elemente mreže
  for (auto const & elem : elements(gridView))
        {
          // Centar elementa
          auto global = elem.geometry().center();
          // Smjesti vrijednost funkcije u centru elementa na odgovarajuće mjesto u vektor.
          //coeff[mapper.index(elem)] = f(global);
          coeff[set.index(elem)] = f(global);
        }
}

#endif /* ELEMENTDATA_HH_ */

**Zadatak**. Potrebno je iscrtati u VTK formatu funkciju i njen interpoland. 
Koristiti `YaspGrid` mrežu i Q1 elemente. Za zadanu funkciju stem:[f] treba iscrtati
f i globalni interpoland Pi_h f. Kako bi se funkcije f i  Pi_h f
mogle razlikovati iscrtavanje treba napraviti na profinjenoj mreži. Na primjer, osnovna mreža
može imati 5 elemenata u 1D ili 5x5 elemenata u 2D, a mreža za iscrtavanje se dobiva sa tri 
profinjenja te mreže (`globalRefine` metoda). 

**Uputa**. Program treba iterirati po svim elementima profinjene mreže (_elementi listovi_) 
i na svakom elementu po svim njegovim  vrhovima. Za element list treba
naći njegov roditelj razine nula (s inicijalne mreže) pomoću metode `father()`. Taj elemenet zovemo 
_makroelement_. Koristeći `geometry` objekt makroelementa možemo vrh elementa lista  na kojem se nalazimo 
preslikati u referentni element i tamo izračunati bazne funkcije. Zatim se  izračunaju koeficijenti 
interpolacijskog operatora i napravi se linearna kombinacija koja predstavlja vrijednost interpolanda u
danom vrhu. Tu vrijednost treba smjestiti u vektor na odgovarajuće mjesto (koristiti `indexSet` profinjene
mreže). 

Pri računu interpolanda funkcije f, funkcija `interpolate` iz klase `LocalInterpolation`
očekuje kao argument funkciju  f\circ g_E, gdje je g_E geometrijsko preslikavanje 
makroelementa na kojem se nalazimo.  Najlakši način da se ostvari kompozicija je da funkcijski objekt
u svom konstruktoru sačuva referencu na makroelement. Na pri
mjer, funkcija koja računa 
sin(2\pi x) bi bila implementirana na sljedeći način:

```c++
template <typename Element, int dim>
struct Function{
    Element const & element;
    Function(Element const & el) : element(el) {}

    void evaluate(Dune::FieldVector<double, dim> const & x_local,
                  Dune::FieldVector<double, 1> & y) const {
        auto x_global = element.geometry().global(x_local);
        evaluate_global(x_global, y);
    }

    void evaluate_global(Dune::FieldVector<double, dim> const & x_global,
                  Dune::FieldVector<double, 1> & y) const {
        y = std::sin(2*M_PI*x_global[0]);
    }
};
```

U petlji po svim elementima novi funkcijski objekt se instancira na svakom elementu 
i predaje mu se pripadni  makroelement. 
Ako je potrebno znati tip makroelementa može se koristiti sljedeći kod:

```c++
using GV = GridType::LeafGridView;
using Element = GV::template Codim<0>::Entity;
```
(tip elementa je isti na svakoj razini). Korištenjem *auto* deklaracije varijabli taj tip uglavnom ne trebamo koristiti 
eksplicitno. 


Za detalje vidjeti stranicu: [bazne_funkcije.html](https://web.math.pmf.unizg.hr/nastava/nrpdj/html/bazne_funkcije.html).

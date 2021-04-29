**Zadatak**. Potrebno je učitati ravninsku mrežu sastavljenu od trokuta
`poluvijenac.msh` (u `gmsh` formatu) i izračunati maksimalno odstupanje od 
ortogonalnosti mreže. 

Reći ćemo da je triangulacija ortogonalna ako svaka dva susjedna trokuta 
dijele čitavu stranicu (triangulacija je konformna) i spojnica njihovih centara
je okomita na stranicu koja ih spaja. 

Nas zanima koliko je odstupanje dane mreže od ortogonalnosti mjereno u stupnjevima.
Za svaku stranicu triangulacije treba izračunati kut _phi_ između spojnice centara 
susjednih elemenata i stranice te odstupanje od ortogonalnosti 
koje definiramo kao 90 - _phi_. Maksimalno
odstupanje mreže od ortogonalnosti je maksimum vrijednosti na svim stranicama. 

Ako je stranica na rubu domene, onda se računa odstupanje od ortogonalnosti između
stranice i spojnice centra elementa i centra stranice. 

Treba napisati program koji učitava mrežu iz .msh datoteke i ispisuje maksimalno odstupanje
od ortogonalnosti mreže. 



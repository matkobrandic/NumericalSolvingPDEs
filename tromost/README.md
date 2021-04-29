*Zadatak*.
Izračunajte moment tromosti 2D tijela zadanog triangulacijom pinch-2D-simplex.msh
oko osi
<a href="https://www.codecogs.com/eqnedit.php?latex=\inline&space;x_1=0" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\inline&space;x_1=0" title="x_1=0" /></a>.
Tijelo ima uniformnu gustoću mase ρ=1kg/cm3. Formula za moment tromosti je

<a href="https://www.codecogs.com/eqnedit.php?latex=I&space;=&space;\int_{\Omega}&space;d(\mathbf{x})^2&space;\rho\,&space;dV" target="_blank"><img src="https://latex.codecogs.com/gif.latex?I&space;=&space;\int_{\Omega}&space;d(\mathbf{x})^2&space;\rho\,&space;dV" title="I = \int_{\Omega} d(\mathbf{x})^2 \rho\, dV" /></a>

gdje je
<a href="https://www.codecogs.com/eqnedit.php?latex=\inline&space;d(\mathbf{x})" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\inline&space;d(\mathbf{x})" title="d(\mathbf{x})" /></a>
 udaljenost točke <a href="https://www.codecogs.com/eqnedit.php?latex=\inline&space;x" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\inline&space;x" title="x" /></a> od osi
<a href="https://www.codecogs.com/eqnedit.php?latex=\inline&space;x_1=0" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\inline&space;x_1=0" title="x_1=0" /></a>.
Koristiti 2D mrežu danu u datoteci `src/pinch-2D-simplex.msh` i `UGGrid`. Pored tromosti treba još izračunati
(radi provjere koda) koordinate težišta tijela. Sve tri vrijednosti će vratiti funkcija `tromost()`
u obliku polja dimenzije 3. Prve dvije vrijednosti su koordinate težišta, a treća je tražena tromost.

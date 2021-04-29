# Konačne diferencije za eliptičku zadaću u 1D

Treba riješiti rubnu zadaću:

<a href="https://www.codecogs.com/eqnedit.php?latex=-u''&space;=&space;f&space;\quad&space;\text{na&space;}&space;(0,L)\\&space;u(0)&space;=&space;g_0,\;&space;u(L)=g_1" target="_blank"><img src="https://latex.codecogs.com/gif.latex?-u''&space;=&space;f&space;\quad&space;\text{na&space;}&space;(0,L)\\&space;u(0)&space;=&space;g_0,\;&space;u(L)=g_1" title="-u'' = f \quad \text{na } (0,L)\\ u(0) = g_0,\; u(L)=g_1" /></a>

 metodom konačnih diferencija s uniformnim prostornim korakom h.
 Diferencijska jednadžba glasi:

<a href="https://www.codecogs.com/eqnedit.php?latex=-(u_{i&plus;1}&space;-&space;2u_i&space;&plus;&space;u_{i-1}&space;)&space;=&space;h^2&space;f_i,\quad&space;i=1,\ldots,n" target="_blank"><img src="https://latex.codecogs.com/gif.latex?-(u_{i&plus;1}&space;-&space;2u_i&space;&plus;&space;u_{i-1}&space;)&space;=&space;h^2&space;f_i,\quad&space;i=1,\ldots,n" title="-(u_{i+1} - 2u_i + u_{i-1} ) = h^2 f_i,\quad i=1,\ldots,n" /></a>

 na mreži <a href="https://www.codecogs.com/eqnedit.php?latex=x_i&space;=&space;i&space;h,&space;h=L/(n&plus;1)" target="_blank"><img src="https://latex.codecogs.com/gif.latex?x_i&space;=&space;i&space;h,&space;h=L/(n&plus;1)" title="x_i = i h, h=L/(n+1)" /></a>. Uočimo da smo jednadžbu diskretizirali samo u unutarnjim točkama domene.
Da bismo uvažili rubne uvjete dodat ćemo dvije trivijalne jednadžbe koje fiksiraju rubni uvjet: 

<a href="https://www.codecogs.com/eqnedit.php?latex=u_0&space;=&space;g_0&space;\text{&space;za&space;}&space;i=0,&space;\quad&space;u_{n&plus;1}&space;=&space;g_1&space;\text{&space;za&space;}&space;i=n&plus;1" target="_blank"><img src="https://latex.codecogs.com/gif.latex?u_0&space;=&space;g_0&space;\text{&space;za&space;}&space;i=0,&space;\quad&space;u_{n&plus;1}&space;=&space;g_1&space;\text{&space;za&space;}&space;i=n&plus;1" title="u_0 = g_0 \text{ za } i=0, \quad u_{n+1} = g_1 \text{ za } i=n+1" /></a>


*Napomena:* Ovakav tretman rubnog uvjeta narušava simetriju matrice.


Rješenje testiramo na funkciji 

<a href="https://www.codecogs.com/eqnedit.php?latex=u&space;=&space;\frac{5L^2}{4\pi^2}\sin(\frac{2\pi&space;x}{L})&space;&plus;&space;(g_1-g_0)\frac{x}{L}&space;&plus;&space;g_0" target="_blank"><img src="https://latex.codecogs.com/gif.latex?u&space;=&space;\frac{5L^2}{4\pi^2}\sin(\frac{2\pi&space;x}{L})&space;&plus;&space;(g_1-g_0)\frac{x}{L}&space;&plus;&space;g_0" title="u = \frac{5L^2}{4\pi^2}\sin(\frac{2\pi x}{L}) + (g_1-g_0)\frac{x}{L} + g_0" /></a>,

koja daje desnu stranu 

<a href="https://www.codecogs.com/eqnedit.php?latex=f&space;=&space;5\sin(\frac{2\pi&space;x}{L})" target="_blank"><img src="https://latex.codecogs.com/gif.latex?f&space;=&space;5\sin(\frac{2\pi&space;x}{L})" title="f = 5\sin(\frac{2\pi x}{L})" /></a>.

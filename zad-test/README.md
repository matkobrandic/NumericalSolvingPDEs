*Zadatak*
U ovom zadatku se konstruira UGGrid čitanjem .msh datoteke.
Potrebno je grid ispisati u VTK formatu.

U zadatku je dovoljno otkomentirati kod na kraju main funkcije.
Zadatak služi za testiranje rada s GitHub repozitorijem.

*Git* je sustav za kontrolu koda (eng. _version control system_). 
Funkcionira tako da za izvorni kod kreira repozitorij u kojem
pamti sve verzije koda koje su spremljene u njega tokom razvoja koda.
 Na taj se način prilikom razvoja koda uvijek
možemo vratiti na prethodnu verziju ukoliko smo unijeli greške 
izmjenama koda. Git omogućava suradnju više programera na istom kodu. 

*GitHub* je servis koji nudi javne i besplatne git repozitorije. Git repozitoriju
koji je smješten na GitHubu može se pristupiti sa bilo kojeg mjesta preko mreže. 
Za rad sa GitHubom svaki korisnik mora kreirati svoj [GitHub _account_](https://github.com/).
Za to je potrebno napraviti _Sign Up_ pri čemu treba dati svoju email adresu i 
odabrati korisničko ime i lozinku.  Zatim treba odabrati besplatne javne
repozitorije jer se privatni repozitoriji u principu plaćaju. 


Naš ciklus rada sastoji se u sljedećem: *GitHub Classroom* (servis GitHub-a za
korištenje repozitorija u nastavi) svakom studentu kreira (besplatan) privatni repozitorij sa 
početnim kodom zadatka. Taj repozitorij treba zatim *klonirati* na lokalni stroj,
napraviti na njemu izmjene i dorade te ga gurnuti na GitHub (u danom roku). 
Nastavnik nakon toga klonira sve studentske repozitorije radi ocjenjivanja zadataka.

Za rad sa GitHubom potrebno je poznavati osnove rada sa Git-om. Za to je 
dovoljno pogledati poglavlja  "Getting Started" i "Git Basics"
u knjizi [Pro Git](https://git-scm.com/book/en/v2).

Prije rada s Git-om potrebno je konfigurirati ga dajući mu svoje ime i prezime te email adresu:

```bash
git config --global user.name "John Doe"
git config --global user.email johndoe@example.com
```

Prvi korak je kloniranje repozitorija na vaš lokalni disk. Treba se pozicionirati 
na mjesto gdje se novi direktorij želi kreirati i otipkate u terminalu

```bash
git clone https://github.com/PMF-MO-NRPDJ/ime_zadatka-username.git ime_zadatka
```

Za svaki zadatak studentu se šalje  poruka putem _Merlin_ sustava. Student je dužan slijediti link 
u poruci nakon čega mu GitHub kreira privatni repozitorij sa zadatkom. Potrebno je ulogirati se u GitHub 
i potražiti adresu koju treba koristiti pri kloniranju na vrhu web-stranice repozitorija.
Svi studentski repozitoriji inaju ime oblika *ime_zadatka-username.git*. Kloniranje treba 
izvršiti u direktorij koji se zove *ime_zadatka* (bez nastavka koji je dan imenom korisnika) 
jer se inače program neće kompilirati. 


Nakon što je repozitorij kloniran potrebno je raditi na zadatku, testirati ga i konačno izmjene 
vratiti na GitHub. 
Kada smo u lokalnoj kopiji zadatka na našem disku izmjenili pojedine datoteke trebamo u 
korijenskom direktoriju zadatka otipkati 

```bash
git status
```

(sve git naredbe počinju sa _git_). To će nam izlistati sve datoteke koje su promijenjene, ali 
čije promjene još nisu registrirane u repozitoriju. Prvi je korak registrirati izmjenjene datoteke 
pomoću naredbe _git add_. Na primjer, ako je _git status_ pokazao samo 

```bash
On branch master
Your branch is up-to-date with 'origin/master'.
Changes not staged for commit:
  (use "git add <file>..." to update what will be committed)
  (use "git checkout -- <file>..." to discard changes in working directory)

        modified:   README.md

Untracked files:
  (use "git add <file>..." to include in what will be committed)

        .README.md.swp

no changes added to commit (use "git add" and/or "git commit -a")
```
onda to znači da je izmijenjena samo datoteku README.md. Datoteka .README.md.swp nije u repozitoriju 
jer je to privremena datoteka koju kreira editor; Nju ignoriramo.

Zatim registriramo naše izmjene pomoću:

```bash
git add README.md
```

(Poslije _git add_ može ići bilo koji broj datoteka.)
Nakon toga _git status_ daje malo drugačiju poruku:

```bash
On branch master
Your branch is up-to-date with 'origin/master'.
Changes to be committed:
  (use "git reset HEAD <file>..." to unstage)

        modified:   README.md

Untracked files:
  (use "git add <file>..." to include in what will be committed)

        .README.md.swp
```

Sad repozitorij *zna* za izmjene u datoteci README.md ali te izmjene još nisu u samom repozitoriju.
Da bismo izmjene ubacijli u repozitorij trebamo otipkati


```bash
git commit -m "Kratko objašnjenje izmjena"
```

sada su izmjene u repozitoriju. Ovaj se ciklus (git add, git commit) može ponoviti
više puta. U _git commit_ naredbi s opcijom -m dajemo tekst kojim opisujemo aše izmjene.
Kada je zadatak gotov tekst izmjene može biti jednostavno "Gotov.".

Zadnji korak je gurnuti izmjene na GitHub. To se radi s naredbom

```bash
git push origin master
```

To je otprilike sve što nam treba za rad sa Git-om i GitHub-om. Više detalja možete naći u knjizi 
  [Pro Git](https://git-scm.com/book/en/v2) i na Internetu. 

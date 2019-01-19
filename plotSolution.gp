n=999
reset
set term gif animate
set output "animate.gif"

set xrange [0:2]
set yrange [-100:500]

set style data lines
do for [i=0:n] {
  plot sprintf('Resultat/sol%i.dat', i) using 1:2 title sprintf("U0 n=%i",i), sprintf('Resultat/sol%i.dat', i) using 1:3 title sprintf("U1 n=%i",i), sprintf('Resultat/sol%i.dat', i) using 1:6 title sprintf("ut n=%i",i),
;}

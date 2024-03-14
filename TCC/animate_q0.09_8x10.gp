reset 
 set term gif animate
set output "animate_q0.09_8x10.gif"
set palette defined (0.25 'black', 0.5 'green',1 'blue')
set cbrange [0.2:1]
set style fill transparent solid 0.9
unset key
set yrange[-12:12]
set xrange[-2:27]

set ylabel "y"
set xlabel "x"
set cblabel "normalized radii"

load 'sequencia.plot'


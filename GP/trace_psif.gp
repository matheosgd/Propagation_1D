#offset1 = 109737.3157383409
#scale   = 0.83000

set term qt font "Times, 15"

set xlabel 'x (grille) [ua]'
set ylabel '||psif(x)||^2'
set title  'paquet d onde à l état final'
set yrange [0:1000000]

plot   '/home/msegaud/stageicp/TD_Matheo/psif' u 1:($8*1 + 186493.14416833327) w l lw 2 t 'psif'
replot '/home/msegaud/stageicp/TD_Matheo/pot' u 2:5 w l lw 2 t 'pot Double Puits'
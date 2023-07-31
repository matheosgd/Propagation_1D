set xlabel 'temps t [ua]'
set ylabel 'Position moyenne [ua]'
set title  'Position moyenne au cours du temps'
#set yrange [0:15000000]

set term qt font "Times, 15"

set xlabel 't [ua]'
set ylabel 'X(t) [ua]'
set title  'Evolution temporelle de la position moyenne'

plot   '/home/msegaud/stageicp/TD_Matheo/prop' u 1:2 w l lw 2 t 'X(t)' #trace Position Moyenne = f(t)
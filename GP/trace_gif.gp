set term gif
set term gif font "Times, 10"
set term gif animate delay 25
set output "output.gif"
scale = 1.5
offset = 186493.14416833327
time = 0
stats 'psif' nooutput

set xrange [-2:2]
set yrange [0:1000000]

set xlabel 'x (grille) [ua]'
set ylabel 'E_p [cm^{-1}]'
set title  'Double Puits; mass=20ua; x0=-1.2ua; sigma=0.4; prop spectral'
set grid

plot   '/home/msegaud/stageicp/TD_Matheo/psif' u 1:($2 *scale + 186493.14416833327) w l lw 2 t '||psi(x,t=0.000ua)||^2', '/home/msegaud/stageicp/TD_Matheo/pot' u 2:5 w l lw 2 t 'pot Double Puits'
plot   '/home/msegaud/stageicp/TD_Matheo/psif' u 1:($3 *scale + 186493.14416833327) w l lw 2 t '||psi(x,t=400ua)||^2'  , '/home/msegaud/stageicp/TD_Matheo/pot' u 2:5 w l lw 2 t 'pot Double Puits'
plot   '/home/msegaud/stageicp/TD_Matheo/psif' u 1:($4 *scale + 186493.14416833327) w l lw 2 t '||psi(x,t=800ua)||^2'  , '/home/msegaud/stageicp/TD_Matheo/pot' u 2:5 w l lw 2 t 'pot Double Puits'
plot   '/home/msegaud/stageicp/TD_Matheo/psif' u 1:($5 *scale + 186493.14416833327) w l lw 2 t '||psi(x,t=1200ua)||^2' , '/home/msegaud/stageicp/TD_Matheo/pot' u 2:5 w l lw 2 t 'pot Double Puits'
plot   '/home/msegaud/stageicp/TD_Matheo/psif' u 1:($6 *scale + 186493.14416833327) w l lw 2 t '||psi(x,t=1600ua)||^2' , '/home/msegaud/stageicp/TD_Matheo/pot' u 2:5 w l lw 2 t 'pot Double Puits'
plot   '/home/msegaud/stageicp/TD_Matheo/psif' u 1:($7 *scale + 186493.14416833327) w l lw 2 t '||psi(x,t=2000ua)||^2' , '/home/msegaud/stageicp/TD_Matheo/pot' u 2:5 w l lw 2 t 'pot Double Puits'
plot   '/home/msegaud/stageicp/TD_Matheo/psif' u 1:($8 *scale + 186493.14416833327) w l lw 2 t '||psi(x,t=2400ua)||^2' , '/home/msegaud/stageicp/TD_Matheo/pot' u 2:5 w l lw 2 t 'pot Double Puits'
plot   '/home/msegaud/stageicp/TD_Matheo/psif' u 1:($9 *scale + 186493.14416833327) w l lw 2 t '||psi(x,t=2800ua)||^2' , '/home/msegaud/stageicp/TD_Matheo/pot' u 2:5 w l lw 2 t 'pot Double Puits'
plot   '/home/msegaud/stageicp/TD_Matheo/psif' u 1:($10*scale + 186493.14416833327) w l lw 2 t '||psi(x,t=3200ua)||^2' , '/home/msegaud/stageicp/TD_Matheo/pot' u 2:5 w l lw 2 t 'pot Double Puits'
plot   '/home/msegaud/stageicp/TD_Matheo/psif' u 1:($11*scale + 186493.14416833327) w l lw 2 t '||psi(x,t=3600ua)||^2' , '/home/msegaud/stageicp/TD_Matheo/pot' u 2:5 w l lw 2 t 'pot Double Puits'
plot   '/home/msegaud/stageicp/TD_Matheo/psif' u 1:($12*scale + 186493.14416833327) w l lw 2 t '||psi(x,t=4000ua)||^2' , '/home/msegaud/stageicp/TD_Matheo/pot' u 2:5 w l lw 2 t 'pot Double Puits'
plot   '/home/msegaud/stageicp/TD_Matheo/psif' u 1:($13*scale + 186493.14416833327) w l lw 2 t '||psi(x,t=4400ua)||^2' , '/home/msegaud/stageicp/TD_Matheo/pot' u 2:5 w l lw 2 t 'pot Double Puits'
plot   '/home/msegaud/stageicp/TD_Matheo/psif' u 1:($14*scale + 186493.14416833327) w l lw 2 t '||psi(x,t=4800ua)||^2' , '/home/msegaud/stageicp/TD_Matheo/pot' u 2:5 w l lw 2 t 'pot Double Puits'
plot   '/home/msegaud/stageicp/TD_Matheo/psif' u 1:($15*scale + 186493.14416833327) w l lw 2 t '||psi(x,t=5200ua)||^2' , '/home/msegaud/stageicp/TD_Matheo/pot' u 2:5 w l lw 2 t 'pot Double Puits'
plot   '/home/msegaud/stageicp/TD_Matheo/psif' u 1:($16*scale + 186493.14416833327) w l lw 2 t '||psi(x,t=5600ua)||^2' , '/home/msegaud/stageicp/TD_Matheo/pot' u 2:5 w l lw 2 t 'pot Double Puits'
plot   '/home/msegaud/stageicp/TD_Matheo/psif' u 1:($17*scale + 186493.14416833327) w l lw 2 t '||psi(x,t=6000ua)||^2' , '/home/msegaud/stageicp/TD_Matheo/pot' u 2:5 w l lw 2 t 'pot Double Puits'

#do for [COL=2:12] { 
#    if (COL=2)  {time = 0}
#    if (COL=3)  {time = 1000}
#    if (COL=4)  {time = 1500}
#    if (COL=5)  {time = 2500}
#    if (COL=6)  {time = 3500}
#    if (COL=7)  {time = 4500}
#    if (COL=8)  {time = 5000}
#    if (COL=9)  {time = 6000}
#    if (COL=10) {time = 6500}
#    if (COL=11) {time = 7000}
#    if (COL=12) {time = 8000}
#    plot     '/home/msegaud/stageicp/TD_Matheo/psif' u 1:(int(COL)*scale + offset) w l lw 2 t sprintf("psi(t = %d)", time), '/home/msegaud/stageicp/TD_Matheo/pot'  u 2:5 w l lw 2 t 'pot Double Puits'
#    }


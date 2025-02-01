#offset1 = 109737.3157383409                    #2047.5923792483 #(morse)  #109737.3157383409 #(harmo) #POUR LES 3 PREMIERS VP
#offset2 = 329211.9479705435                    #6012.7138648988           #329211.9479705435
#offset3 = 548686.5970521740                    #9804.4852272820           #548686.5970521740
#scale   = 0.83000                              #0.00625                   #0.83000

#set xlabel 'x (grille) [ua]'
#set ylabel 'E_p [cm^{-1}]'
#set title  'Paquet d onde à tf = 500.0 ua sur pot. Morse'# et trois premiers etats propres associes'

set term qt font "Times, 15"

#plot '/home/msegaud/stageicp/TD_Matheo/pot'  u 2:5 w l lw 2 t 'pot Harmonique'

#TRACE LE PAQUET D'ONDE INITIAL + DEUX TEMPS DIFFERENTS
#scale = 22.5
#set title 'Paquet d ondes à trois t différents (pot. Harmonique)'# et trois premiers etats propres associes'
#set yrange [0:60000]
#replot   '/home/msegaud/stageicp/TD_Matheo/mat' u 1:(abs($2)*scale+797424.62009915092) w l lw 2 t 'Psi(t = t0)' 
#replot '/home/msegaud/stageicp/TD_Matheo/mat' u 1:(abs($3)*scale+797424.62009915092) w l lw 2 t 'Psi(t = t1)' 
#replot '/home/msegaud/stageicp/TD_Matheo/mat' u 1:(abs($4)*scale+797424.62009915092) w l lw 2 t 'Psi(t = tf)' 

#TRACE TROIS PREMIERS ETATS PROPRES ET POTENTIEL
#replot   '/home/msegaud/stageicp/TD_Matheo/mat' u 1:($2*scale+offset1) w l lw 2 t 'Eigvec_1' 
#replot '/home/msegaud/stageicp/TD_Matheo/mat' u 1:($3*scale+offset2) w l lw 2 t 'Eigvec_2'
#replot '/home/msegaud/stageicp/TD_Matheo/mat' u 1:($4*scale+offset3) w l lw 2 t 'Eigvec_3'

#TRACE SPECIAL HARMONIQUE /!\ mettre nb = 15; nq = 30
#set title  'Paquet d ondes à tf = 500.0 ua sur pot. Harmonique'
#scale = 22.5
#set yrange [0:4000000]

#replot   '/home/msegaud/stageicp/TD_Matheo/psif' u 1:(abs($2)*scale+797424.62009915092) w l lw 2 t 'spectral' 
#replot '/home/msegaud/stageicp/TD_Matheo/psif' u 1:(abs($3)*scale+797424.62009915092) w l lw 2 t 'euler' 
#replot '/home/msegaud/stageicp/TD_Matheo/psif' u 1:(abs($4)*scale+797424.62009915092) w l lw 2 t 'taylor_4' 
#replot '/home/msegaud/stageicp/TD_Matheo/psif' u 1:(abs($5)*scale+797424.62009915092) w l lw 2 t 'SIL' 
#replot '/home/msegaud/stageicp/TD_Matheo/pot'  u 2:5 w l lw 2 t 'pot Harmonique'

#TRACE SPECIAL MORSE
#set title  'Paquet d ondes à tf = 500.0 ua sur pot. Morse'
#scale = 0.05
#set yrange [0:65000]
#replot   '/home/msegaud/stageicp/TD_Matheo/psif' u 1:(abs($2)*scale+14881.271632364544) w l lw 2 t 'spectral' 
#replot '/home/msegaud/stageicp/TD_Matheo/psif' u 1:(abs($3)*scale+14881.271632364544) w l lw 2 t 'euler' 
#replot '/home/msegaud/stageicp/TD_Matheo/psif' u 1:(abs($4)*scale+14881.271632364544) w l lw 2 t 'taylor_4' 
#replot '/home/msegaud/stageicp/TD_Matheo/psif' u 1:(abs($5)*scale+14881.271632364544) w l lw 2 t 'SIL' 
#replot '/home/msegaud/stageicp/TD_Matheo/pot'  u 2:5 w l lw 2 t 'pot Morse'

#TRACE ANIMATION
scale = 1.5
set xlabel 'x (grille) [ua]'
set ylabel '||psif(x)||^2'
set title  'paquet d onde à l état final'
set yrange [0:1000000]

plot     '/home/msegaud/stageicp/TD_Matheo/psif' u 1:($2 *scale + 186493.14416833327) w l lw 2 t 'psi0'
replot   '/home/msegaud/stageicp/TD_Matheo/psif' u 1:($3 *scale + 186493.14416833327) w l lw 2 t 'psif1'
replot   '/home/msegaud/stageicp/TD_Matheo/psif' u 1:($4 *scale + 186493.14416833327) w l lw 2 t 'psif2'
replot   '/home/msegaud/stageicp/TD_Matheo/psif' u 1:($5 *scale + 186493.14416833327) w l lw 2 t 'psif3'
replot   '/home/msegaud/stageicp/TD_Matheo/psif' u 1:($6 *scale + 186493.14416833327) w l lw 2 t 'psif4'
replot   '/home/msegaud/stageicp/TD_Matheo/psif' u 1:($7 *scale + 186493.14416833327) w l lw 2 t 'psif5'
replot   '/home/msegaud/stageicp/TD_Matheo/psif' u 1:($8 *scale + 186493.14416833327) w l lw 2 t 'psif6'
replot   '/home/msegaud/stageicp/TD_Matheo/psif' u 1:($9 *scale + 186493.14416833327) w l lw 2 t 'psif7'
replot   '/home/msegaud/stageicp/TD_Matheo/psif' u 1:($10*scale + 186493.14416833327) w l lw 2 t 'psif8'
replot   '/home/msegaud/stageicp/TD_Matheo/psif' u 1:($11*scale + 186493.14416833327) w l lw 2 t 'psif9'
replot   '/home/msegaud/stageicp/TD_Matheo/psif' u 1:($12*scale + 186493.14416833327) w l lw 2 t 'psif10'
replot   '/home/msegaud/stageicp/TD_Matheo/psif' u 1:($13*scale + 186493.14416833327) w l lw 2 t 'psif11'
replot   '/home/msegaud/stageicp/TD_Matheo/psif' u 1:($14*scale + 186493.14416833327) w l lw 2 t 'psif12'
replot   '/home/msegaud/stageicp/TD_Matheo/psif' u 1:($15*scale + 186493.14416833327) w l lw 2 t 'psif13'
replot   '/home/msegaud/stageicp/TD_Matheo/psif' u 1:($16*scale + 186493.14416833327) w l lw 2 t 'psif14'
replot   '/home/msegaud/stageicp/TD_Matheo/psif' u 1:($17*scale + 186493.14416833327) w l lw 2 t 'psif15'


replot   '/home/msegaud/stageicp/TD_Matheo/pot'  u 2:5 w l lw 2 t 'pot Double Puits'

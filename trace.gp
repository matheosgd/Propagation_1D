#offset1 = 109737.3157383409                    #2047.5923792483 #(morse)  #109737.3157383409 #(harmo) #POUR LES 3 PREMIERS VP
#offset2 = 329211.9479705435                    #6012.7138648988           #329211.9479705435
#offset3 = 548686.5970521740                    #9804.4852272820           #548686.5970521740
#scale   = 0.83000                              #0.00625                   #0.83000

set xlabel 'x (grille) [ua]'
set ylabel 'E_p [cm^{-1}]'
set title  'Paquet d onde à tf = 500.0 ua sur pot. Morse'# et trois premiers etats propres associes'
set yrange [0:65000]

set term qt font "Times, 15"

#TRACE LE PAQUET D'ONDE INITIAL + DEUX TEMPS DIFFERENTS
#scale = 22.5
#set title 'Paquet d ondes à trois t différents (pot. Harmonique)'# et trois premiers etats propres associes'
#set yrange [0:60000]
#plot   '/home/msegaud/stageicp/TD_Matheo/mat' u 1:(abs($2)*scale+797424.62009915092) w l lw 2 t 'Psi(t = t0)' 
#replot '/home/msegaud/stageicp/TD_Matheo/mat' u 1:(abs($3)*scale+797424.62009915092) w l lw 2 t 'Psi(t = t1)' 
#replot '/home/msegaud/stageicp/TD_Matheo/mat' u 1:(abs($4)*scale+797424.62009915092) w l lw 2 t 'Psi(t = tf)' 

#TRACE TROIS PREMIERS ETATS PROPRES ET POTENTIEL
#plot   '/home/msegaud/stageicp/TD_Matheo/mat' u 1:($2*scale+offset1) w l lw 2 t 'Eigvec_1' 
#replot '/home/msegaud/stageicp/TD_Matheo/mat' u 1:($3*scale+offset2) w l lw 2 t 'Eigvec_2'
#replot '/home/msegaud/stageicp/TD_Matheo/mat' u 1:($4*scale+offset3) w l lw 2 t 'Eigvec_3'

#TRACE SPECIAL HARMONIQUE /!\ mettre nb = 15; nq = 30
#set title  'Paquet d onde à tf = 500.0 ua sur pot. Harmonique'
#scale = 22.5
#plot   '/home/msegaud/stageicp/TD_Matheo/psif' u 1:(abs($2)*scale+797424.62009915092) w l lw 2 t 'spectral' 
#replot '/home/msegaud/stageicp/TD_Matheo/psif' u 1:(abs($3)*scale+797424.62009915092) w l lw 2 t 'euler' 
#replot '/home/msegaud/stageicp/TD_Matheo/psif' u 1:(abs($4)*scale+797424.62009915092) w l lw 2 t 'taylor_4' 
#replot '/home/msegaud/stageicp/TD_Matheo/psif' u 1:(abs($5)*scale+797424.62009915092) w l lw 2 t 'SIL' 
#replot '/home/msegaud/stageicp/TD_Matheo/pot'  u 2:5 w l lw 2 t 'pot Harmonique'

#TRACE SPECIAL MORSE
#set title  'Paquet d onde à tf = 500.0 ua sur pot. Morse'
#scale = 0.05
#plot   '/home/msegaud/stageicp/TD_Matheo/psif' u 1:(abs($2)*scale+14881.271632364544) w l lw 2 t 'spectral' 
#replot '/home/msegaud/stageicp/TD_Matheo/psif' u 1:(abs($3)*scale+14881.271632364544) w l lw 2 t 'euler' 
#replot '/home/msegaud/stageicp/TD_Matheo/psif' u 1:(abs($4)*scale+14881.271632364544) w l lw 2 t 'taylor_4' 
#replot '/home/msegaud/stageicp/TD_Matheo/psif' u 1:(abs($5)*scale+14881.271632364544) w l lw 2 t 'SIL' 
#replot '/home/msegaud/stageicp/TD_Matheo/pot'  u 2:5 w l lw 2 t 'pot Morse'

plot '/home/msegaud/stageicp/TD_Matheo/pot'  u 2:5 w l lw 2 t 'pot Harmonique'

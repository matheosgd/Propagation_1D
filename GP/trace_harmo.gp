set term qt font "Times, 15"

set xlabel 'x (grille) [ua]'
set ylabel 'E_p [cm^{-1}]'
set title  'Pot. Harmonique'                  #et trois 1^{ers} états propres associés'
set yrange [0:1000000]

#TRACE POT ET TROIS PREMIERS ETATS PROPRES
#offset1 = 109737.3157383409                  #109737.3157383409 #(harmo) #POUR LES 3 PREMIERS VP
#offset2 = 329211.9479705435                  #329211.9479705435
#offset3 = 548686.5970521740                  #548686.5970521740
#scale   = 0.83000                            #0.83000
#plot   '/home/msegaud/stageicp/TD_Matheo/mat' u 1:($2*scale+offset1) w l lw 2 t 'Eigvec_1' 
#replot '/home/msegaud/stageicp/TD_Matheo/mat' u 1:($3*scale+offset2) w l lw 2 t 'Eigvec_2'
#replot '/home/msegaud/stageicp/TD_Matheo/mat' u 1:($4*scale+offset3) w l lw 2 t 'Eigvec_3'

plot '/home/msegaud/stageicp/TD_Matheo/pot' u 2:5 w l lw 2 t 'pot Harmonique'


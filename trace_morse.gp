set term qt font "Times, 15"

set xlabel 'x (grille) [ua]'
set ylabel 'E_p [cm^{-1}]'
set title  'Pot. de Morse et trois'           #1^{ers} états propres associés'
set yrange [0:12000]

#TRACE POT ET TROIS PREMIERS ETATS PROPRES
#offset1 = 2047.5923792483                    #2047.5923792483 #(morse)   #POUR LES 3 PREMIERS VP
#offset2 = 6012.7138648988                    #6012.7138648988           
#offset3 = 9804.4852272820                    #9804.4852272820           
#scale   = 0.00625                            #0.00625           
#plot   '/home/msegaud/stageicp/TD_Matheo/mat' u 1:($2*scale+offset1) w l lw 2 t 'Eigvec_1' 
#replot '/home/msegaud/stageicp/TD_Matheo/mat' u 1:($3*scale+offset2) w l lw 2 t 'Eigvec_2'
#replot '/home/msegaud/stageicp/TD_Matheo/mat' u 1:($4*scale+offset3) w l lw 2 t 'Eigvec_3'

plot '/home/msegaud/stageicp/TD_Matheo/pot' u 2:5 w l lw 2 t 'pot Morse'
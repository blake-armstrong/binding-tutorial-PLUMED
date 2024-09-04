set terminal pngcairo size 1200,800
set key noenhanced
set output 'fes_radiuses.png'

s = system("ls radius_0.?/aligned_fes.dat") 
set yrange [-32:20]
set xlabel "z-distance (nm)" font "Helvetica,18"
set ylabel "PMF (kJ/mol)" font "Helvetica,18"
set tics font "Helvetica,18"
set key font "Helvetica,18"
set key right bottom
set border lw 3
set xtics scale 2
set ytics scale 2
p "radius_0.0/aligned_fes.dat" u 1:2 w l title "radius = 0.0 nm" lw 2, \
  "radius_0.1/aligned_fes.dat" u 1:2 w l title "radius = 0.1 nm" lw 2, \
  "radius_0.2/aligned_fes.dat" u 1:2 w l title "radius = 0.2 nm" lw 2, \
  "radius_0.3/aligned_fes.dat" u 1:2 w l title "radius = 0.3 nm" lw 2, \
  "radius_0.4/aligned_fes.dat" u 1:2 w l title "radius = 0.4 nm" lw 2, \
  "radius_0.5/aligned_fes.dat" u 1:2 w l title "radius = 0.5 nm" lw 2, \
  -30*exp(-((x)**2/(2*(0.1**2))))+15*exp(-((x-0.3)**2/(2*(0.1**2))))+14473.0*0.5*((x < -0.125)?x+0.125:0)**2 title "Input PES" lw 4 dt "- "

set samples 100
set isosamples 100
set xyplane 0

set xlabel "STRIDES" font ",10"
set ylabel "SIZE (KB, log_2)" font ",10"
set zlabel "MBs/sec" offset -3, 0 font ",10" rotate parallel

set tics font ",8"

set pm3d
splot [:15][17:4] "output.txt" u 1:(log($2)/log(2)):3 matrix nonuniform with lines lc 0 notitle


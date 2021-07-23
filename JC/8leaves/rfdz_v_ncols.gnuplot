set key bottom right

se(n,p) = sqrt(p*(1.0-p)/n)

set xlabel '# columns'
set ylabel 'fraction reconstructed with correct topology'

set style data errorbars

N = 1000
plot \
	'bl0.07_rfds' using 6:10:(se(N,$10)) ps 0.8 t'bl=0.07', \
	'bl0.1_rfds' using 6:10:(se(N,$10)) ps 0.8 t'bl=0.1', \
	'bl0.14_rfds' using 6:10:(se(N,$10)) ps 0.8 t'bl=0.14', \
	'bl0.2_rfds' using 6:10:(se(N,$10)) ps 0.8 t'bl=0.2', \
	'bl0.28_rfds' using 6:10:(se(N,$10)) ps 0.8 t'bl=0.28', \
        'bl0.4_rfds' using 6:10:(se(N,$10)) ps 0.8 t'bl=0.4', \
        'bl0.56_rfds' using 6:10:(se(N,$10)) ps 0.8 t'bl=0.56', \
	0.9 ls 0.5

pause -1

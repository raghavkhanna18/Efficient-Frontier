#Overview

Basic script to calculate the efficient frontier between 2 portfolios with 4 equities of different weights.

Currently, the Variance, Covariance, and Expected Return values are guesstimates. 

#TODO
Calculate Variance, Covariance, and Expected Return values from an input data set.

#Graphing Results using gnuplot
The following commands can be used to plot the results on a graph using gnuplot in terminal:

Assuming the enveloped data is saved to `./envelope.dat` and the efficient_frontier calculations are saved to `./efficient_frontier.dat`

`set key top right`

`set key box`

`set xlabel "Volatility"`

`set ylabel "Expected Return"`

`plot './envelope.dat' using 2:3 with linespoints title "Envelope", './efficient_frontier.dat' using 1:2  with points pointtype 5 title "Efficient Frontier"`
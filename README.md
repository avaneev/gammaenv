# gammaenv #
## Introduction ##

"gammaenv" produces smoothed-out S-curve envelope signal with the specified
attack and release characteristics. The attack and release times can be
further adjusted in real-time. Delay parameter is also specified as the
percentage of the total time.

The S-curve produced by this envelope algorithm closely resembles a
sine-wave signal slightly augmented via the tanh() function. Such
augmentation makes the shape slightly steeper and in the end allows the
algorithm to follow it closer. The name "gammaenv" relates to this
algorithm's version.

The algorithm's topology is based on 5 sets of "leaky integrators" (the
simplest form of 1st order low-pass filters). Each set (except the 5th) use
4 low-pass filters in series. Outputs of all sets are then simply
summed/subtracted to produce the final result. The topology is numerically
stable for any valid input signal, but may produce envelope overshoots
depending on the input signal.

Up to 25% of total attack (or release) time can be allocated (via Delay
parameters) to the delay stage. The delay is implemented by the same
topology: this has the benefit of not requiring additional memory
buffering. For example, it is possible to get up to 250 ms delay on
a 1-second envelope release stage without buffering.

The processSymm() function provides the "perfect" implementation of the
algorithm, but it is limited to the same attack and release times. A more
universal process() function can work with any attack and release times,
but it is about 2 times less efficient and the actual attack stage's
envelope can range from the "designed" U to the undesired sharp V shape.
Unfortunately, the author was unable to find an approach that could be
similar to the processSymm() function while providing differing attack and
release times (the best approach found so far lengthens the delay stage
unpredictably).

Compile and run the testtable.cpp to produce a tab-delimited table of curves
at various delay values. Also use this source as an example of envelope
generator setup and use procedure.

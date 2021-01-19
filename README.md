
# natPsoho

This is an extension of the particle swarm optimization structure
learning algorithm for higher-order dynamic Bayesian networks from
Santos and Maciel (<https://doi.org/10.1109/BRC.2014.6880957>) to an
integer number encoding that supports multiple orders at the same time.
The original algorithm can be found in both the ‘dbnR’ package
(<https://github.com/dkesada/dbnR>) and in its own repository
(<https://github.com/dkesada/PSOHO>). I couldn’t fork my own repository
with the implementation for the original algorithm, so I imported it
into a new one. I will use it as a template and extend it to the integer
number encoding.

The main objective of this new encoding is to be able to learn the
appropriate Markovian order of the network and the structure at the same
time instead of having to fix the order beforehand. In this extension,
you can fix a maximum order and the algorithm will search structures up
to that point. The new encoding offers integer number vectors whose size
only depend on the number of variables in t0, not on the order of the
network. It also offers a way to favour smaller orders by fixing a
geometric distribution that defines how often higher order arcs are
added.

The algorithm is now finished and tested. I represents a big improvement
over the original binary version and it scales much better for higher
orders. Some features are still work in progress, but for the most part
it’s done. The algorithm will be summarized for a conference paper and
afterwards it will be added to the ‘dbnR’ package.

## Some conclusions and remarks

### Similarities with the binary case

The new encoding is similar to the binary case in that the integer
vectors can be translated into binary ones with some maximum order and
it could be argued that both approaches are equivalent. However, the
integer representation of positions and velocities offers the advantage
that there is really no need to fix any order or any maximum order. The
maximum order is useful to perform the search in a promising space
rather than exactly that of the given Markovian order. Also, the
addition of different probabilities to add arcs depending on the order
prevents the network from freely searching in a bigger search space in
favor of graphs with more populated lower orders.

The performance and efficiency are also heavily improved in the new
algorithm. This is mainly due to three reasons:

  - Encoding the structure of the network in vectors whose length is
    quadratic on the number of nodes in *only* t\_0 means that no matter
    the order of the desired network, the size of the particles will
    remain constant. This helps the algorithm greatly when dealing with
    higher orders, where as in the original algorithm the causal lists
    get bigger with each order increase.

  - Performing bitwise operations over an integer is much faster than
    iterating through lists and operating each value independently.

  - Both of my implementations use ‘Rcpp’ and switch to C++ for
    operating over the particles One implementation uses named
    Rcpp::List and the other uses only Rcpp::NumericVector, but when
    implementing the new algorithm as a fork of the old one I saw a lot
    of old parts that could be improved to be more efficient. Those
    inefficiencies are still present in the old algorithm, and so the
    new one is slightly improved in some sections.

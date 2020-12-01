
# intPSOHO

This is an extension of the particle swarm optimization structure
learning algorithm for higher-order dynamic Bayesian networks from
Santos and Maciel (<https://doi.org/10.1109/BRC.2014.6880957>) to
an integer number encoding that supports multiple orders at the
same time. The original algorithm can be found in both the _dbnR_ 
package and in (<https://github.com/dkesada/PSOHO>). I couldn't fork 
my own repository with the implementation for the original algorithm, 
so I imported it into a new one. I will use it as a template and 
extend it to the integer number encoding.

The main objective of this new encoding is to be able to learn the
appropiate Markovian order of the network and the structure at the
same time instead of having to fix the order beforehand. In this
extension, you can fix a maximum order and the algorithm will search
structures up to that point. The new encoding offers integer number
vectors whose size only depend on the number of variables in t0, not
on the order of the network. It also offers a way to favor smaller 
orders by fixing a geometric distribution that defines how often
higher order arcs are added.

The final objectives are to add this new extension to the _dbnR_ 
package and to create a conference article describing it.

## Some conclusions and remarks

### Similarities with the binary case

The new encoding is similar to the binary case in that the integer
vectors can be translated into binary ones with some maximum order
and it could be argued that both approaches are equivalent. 
However, the integer representation of positions and velocities
offers the advantage that there is really no need to fix any order
or any maximum order. The maximum order is useful to perform the 
search in a promising espace rather than exactly that of the 
given Markovian order. Also, the addition of different 
probabilities to add arcs depending on the order prevents the 
network from freely searching in a bigger search space in favor
of graphs with more populated lower orders.

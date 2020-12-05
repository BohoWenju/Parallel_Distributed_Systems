# Parallel/Distributed Systems
Calculation of triangles in an undirected weightless graph implemented both in serial and parallel with two separate algorithms(in c)

  Aikaterini Prokou/Eleftherios Mourelatos


  Calculates the difference in time
  to calculate the number of triangles in a
  given undirected graph.The graphs were
  taken from Matrix market:
  
  (See :https://math.nist.gov/MatrixMarket/)
  and have been formatted according to the CSC Format
  
  V3 : indicates a simple algorithm to calculate the number of triangles in an undirected weightless graph
  
  V4 : indicates a more complex algorithm to calculate the number of triangles
  in an undirected weightless graph
  based on the formula C=(A.*(A*A))*e/2 where : C is a vector with the number of triangles each node is connected to ,
  A is the matrix as taken from matrix market and e is a vector filled with 1's.Number of triangles is then
  the sum of vector C divided by 3.



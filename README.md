# Fast and accurate evaluation of close-to-touching rigit body interations
### Flatiron summer 2023 internship

We consider the flow of dense suspensions of rigid bodies in a Stokesian fluid. Such flows are difficult to compute numerically due to several close-to-touching interactions. Such interactions require a large number of unknowns to resolve sharply peaked surface forces/densities. The number of GMRES iterations required to solve the discretized PDE is also usually large. The time evolution of such flows is very stiff and requires taking extremely small time steps. A common way of dealing with these difficulties is to introduce a repulsion force between the particles, to prevent them from getting too close. However, this additional repulsion force is non-physical and may alter the results of these simulations.

For suspensions of identical discs in 2D, we present a fast and accurate method that tries to mitigate some of these challenges without introducing artificial forces. Our method solves the issue of large number of unknowns in the PDE discretization through precomputation, compression and interpolation of the close-to-touching part of the interaction operator. Our method also helps to improve the convergence rate of the GMRES solves.





#### Timeline

  - Week 1 (29 May - 2 Jun): read about the problem. Familiarize with tools needed.
  - Week 2 (5 -9 Jun): Read more about boundary integral methods, finish onboarding.
  - Week 3 (12 - 16 Jun): Understand chunkie, Green's identity test.
  - Week 4 (19 - 23 Jun): Capacitance and Elastance problems

# MPI
2D heat transfer parallel solver through master-slave scheme.
Solving 2D heat transfer for a rectangular domain. The domin is divided into many subdomains for parallel programming. 
The Master divides the domain and assigns the chunks of the domains to the Slaves. 
The Master sends information regarding the boundary conditions to each Slave.
The Slaves solve PDEs based on their bundary conditions.
The Slaves send back their results to the Master.
The Master combine the information from each subdomain solved by Slaves.
The Master updates the boundary conditions of each subdomain and sends it to the Slaves to solve the new time step.

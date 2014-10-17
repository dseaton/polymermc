polymermc
=========

Monte Carlo simulations of coarse-grained polymers: Metropolis, Wang-Landau for Flexible and Semiflexible chains

A long overdue push to git, polymermc contains source code used to simulate both flexible and semiflexible homopolymer chains. A single chain can be simulated using the Metropolis algorithm or the Wang-Landau algorithm. Published results include:

* Seaton, D. T., Schnabel, S., Landau, D. P., & Bachmann, M. (2013). From Flexible to Stiff: Systematic Analysis of Structural Phases for Single Semiflexible Polymers. Physical review letters, 110(2), 028103.

* Seaton, D. T., Wüst, T., & Landau, D. P. (2010). Collapse transitions in a flexible homopolymer chain: Application of the Wang-Landau algorithm. Physical Review E, 81(1), 011802.

The two papers above also represent the 1D fully flexible source code (1D_source_code) and the 2D semiflexible simulation along the energy and flexibility axes. Other papers have also utilized previous versions of this code. The most detailed information on utilizing this code can be found in my dissertation:

http://athenaeum.libs.uga.edu/handle/10724/27001

Users are welcome to use the code, and better yet, encouraged to improve the code. This code base was started in 2004 when I was teaching myself to program, and the last major additions were made in 2010. In my freetime, I'll be working to make this code more modular. If users are interested in results, I'll be glad to load pieces of the data (note, the 2D code led to data on the order of a TB).

![Flexible Behavior](https://github.com/dseaton/polymermc/blob/master/images/flexible_behavior.pdf)

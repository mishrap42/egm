# egm
A modular implementation of the method of endogenous gridpoints within the context of economics

# Project Description

The method of endogenous gridpoints (EGM) is an algorithm for solving dynamic programming problems popularized largely in 2005 in a National Bureau of Economics Research paper by Christopher D. Carroll (see more here: http://www.nber.org/papers/t0309). While it is largely used in economics, as a methodology it has quite a bit of potential outside of this core field. Given a Bellman Equation to describe a dynamic programming scenario, the core of assumptions of EGM are fairly ubiquitous enough to make EGM a viable solution methodology. Seeing as dynamic programming has been used in physics, dynamical systems, mechanical systems analysis, sociological and epidemiological models, as just a sampling of the many fields, the generalization of this code can provide a useful tool for solving a whole range of problems.


# Implementation

This project takes the classic EGM code of Carroll, adds the "generalized" EGM described by later economists like Barillas, and presents a commented modular object-oriented version of the code. The code is in Matlab, and the only portion of the code which would need to be changed are the definitions of marginal utility (which in other cases would be the first derivative of a constraint function) and its inverse, the function `marg_val_end_per` which contains the iterator calculating the Bellman Equation, and the value of desired parameters. The code as is also incorporates an endogenous labor component which means that the code demonstrates a multivariate output function rather than a univariate one.

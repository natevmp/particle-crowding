# particle-crowding
Solid disk particles undergoing elastic collisions and growth

## File structure
All files related to the simulation code are found in the src folder. The simulation itself is packaged in the BParts module, which is spread accross the following files:

bmparticles.jl

--> main file; defines objects "Cell", "Bounds", and "Arena", and wraps other files in module

arenabuilder.jl

--> contains methods to create arenas and cells

dynamics.jl

--> contains methods to evolve arenas in time

positions.jl

--> contains methods relating to cell positions

movement.jl

--> contains methods relating to cell movement

collisions.jl

--> contains methods relating to cell collisions

growth.jl

--> contains methods related to population growth

scientist.jl

--> contains methods related to analysing arenas, as well as the the function "randArenaEvolve", which is a premade recipe for running a full simulation.

util.jl

--> contains general utility functions called by other methods

Furthermore, the src folder contains the file bmtheory.jl, which defines the module Theorist, containing functions that do not interact with the simulations or its output directly, but are rather related to the analytical theory (our Langevin equation). It contains the theoretical predictions for statistical quantities such as the mean free path, the time dependent friction coefficient, etc. Furthermore it contains simulations of the Langevin equation and functions related to analysing them.

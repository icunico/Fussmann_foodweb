
## What is fussmann_ml ?

The present model is based on the Fussmann food web (https://onlinelibrary.wiley.com/doi/pdfdirect/10.1046/j.1461-0248.2002.00329.x), 
with extensions to include the microbialloop and mass-balance constraints.
The microbial loop is a key process at the lower level of the trophic network, where bacteria facilitate the recycling of detritus
from higher trophic levels into nutrients, which are subsequently available for the growth of primary producers.
The fussmann_ml code computes the system’s dynamic behavior (species extinction, steady state, periodic oscillations, or
chaos) for all user-defined parameter combinations.

The model implementation requires the following Python packages: ordpy, numpy, scipy.integrate, matplotlib.pyplot




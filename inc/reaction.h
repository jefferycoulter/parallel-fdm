#ifndef REACTION_INCL
#define REACTION_INCL

#include <stdio.h>

/**
 * @brief different types of reactions used to determine an appropriate PDE for
 * a given chemical species
 */
enum Reactions{ReactionDiffusion, Diffusion, Reaction};

/**
 * @brief struct containing information relevant to a chemical species
 */
typedef struct
{
    int type; // type of PDE that describes the species (i.e. reaction-diffusion, diffusion, etc.)
    float ic; // initial conditions
    int bc_type; // type of boundary conditions (either Dirichlet or VonNeumann)
    float bc; // boundary conditions
    float *consts; // constants. diffusion is first, followed by rate constants
    char *name; // name of species
    float *u_now; // current result, used to compute next fdm iteration
    float *u_next; // stores new data after fdm iteration
} Species;

/**
 * @brief struct defining a chemical reaction network
 */
typedef struct
{
    int n_species; // number of species
    Species *species; // array of chemical species in the reaction network
} ReactionNetwork;

void CreateReactionNetwork(ReactionNetwork *network, FILE *fp);

void ParseReaction(FILE *fp);

#endif // REACTION_INCL
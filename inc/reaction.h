#ifndef REACTION_INCL
#define REACTION_INCL

#include <stdio.h>

/**
 * @brief different types of reactions used to determine an appropriate PDE for
 * a given chemical species
 * @param RD reaction-diffusion
 * @param D diffusion
 * @param R reaction
 */
enum ReactionTypes{RD, D, R};

/**
 * @brief struct containing information relevant to a chemical species
 */
typedef struct
{
    char name[24]; // name of species (arbitrary since size, assumes species have short names, i.e. Ca, Sr, IP3, etc.
    float *u_now; // current result, used to compute next fdm iteration
    float *u_next; // stores new data after fdm iteration
} Species;

/**
 * @brief struct defining a chemical reaction
 */
typedef struct
{
    int type; // type of PDE that describes the species (i.e. reaction-diffusion, diffusion, etc.)
    float ic; // initial condition value
    int bc_type; // type of boundary conditions (either Dirichlet or VonNeumann)
    float bc; // boundary condition value

    float *consts; // constants. diffusion is first, followed by rate constants

    int n_species; // number of species in this reaction
    Species *species; // array of chemical species in this reaction
} Reaction;

/**
 * @brief struct defining a chemical reaction network
 */
typedef struct
{
    int n_reactions; // number of reactions in the network
    Reaction *reactions; // array of reactions.
} ReactionNetwork;

void CreateReactionNetwork(ReactionNetwork *network, FILE *fp);

FILE *ParseReaction(ReactionNetwork *network, FILE *fp);

#endif // REACTION_INCL
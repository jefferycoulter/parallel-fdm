#include "reaction.h"

#include <stdlib.h>

void CreateReactionNetwork(ReactionNetwork *network, FILE *fp)
{
    // pass in a ReactionNetwork pointer
    // allocate memory for ReactionNetwork
    network = malloc(sizeof(ReactionNetwork));
    network->reactions = (float*)malloc(sizeof(Species) * (*network).n_reactions);
    while (!feof(fp))
    {
        FILE *fp_temp;
        fp_temp = ParseReaction(network, fp);
        fp = fp_temp;
        // allocate solution arrays
        // determine the type of equation to be used to model the species
    }
}

FILE *ParseReaction(ReactionNetwork *network, FILE *fp)
{
    // if species isn't already in network->reactions->species then add it 
    for (int r = 0; r < network->n_reactions; r++)
    {
        // read line of file

        // break after each colon

        // after first colon, save each term before next colon as a species in dictionary 
        // if it doesn't already exist

        // save reaction rate
    }
}
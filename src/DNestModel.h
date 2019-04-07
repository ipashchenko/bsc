#ifndef BSC__DNESTMODEL_H_
#define BSC__DNESTMODEL_H_

#include "SkyModel.h"
#include "Gains.h"
#include "RNG.h"


class DNestModel {
    public:

        DNestModel();

        // Generate the point from the prior
        void from_prior(DNest4::RNG& rng);

        // Metropolis-Hastings proposals
        double perturb(DNest4::RNG& rng);

        // Likelihood function
        double log_likelihood() const;

        // Print to stream
        void print(std::ostream& out) const;

        // Return string with column information
        std::string description() const;

    private:

        SkyModel* sky_model;
        Gains* gains;

};

#endif //BSC__DNESTMODEL_H_

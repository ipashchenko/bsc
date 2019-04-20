#ifndef BSC__DNESTMODEL_H_
#define BSC__DNESTMODEL_H_

#include "SkyModel.h"
#include "Gains.h"
#include "RNG.h"


class DNestModel {
    public:

        DNestModel();
        ~DNestModel();
        DNestModel(const DNestModel& other);
        DNestModel& operator=(const DNestModel& other);

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

        void set_x_skymodel(double x);

    private:

        SkyModel* sky_model{};
        Gains* gains{};
        // Prediction of SkyModel + Gains
        std::valarray<double> mu_real_full;
        std::valarray<double> mu_imag_full;
        // This runs ``ft`` method of SkyModel with (u, v) from Data and updates SkyModel predictions
        void calculate_sky_mu();
        // This calculates full model (SkyModel with gains) and updates ``mu_real/imag_full``
        void calculate_mu();

};

#endif //BSC__DNESTMODEL_H_

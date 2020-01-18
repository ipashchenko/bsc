#ifndef BSC_MODEL_H
#define BSC_MODEL_H

#include <vector>
#include "Component.h"
#include "RNG.h"


class SkyModel {
    public:
        explicit SkyModel(int n_comp);
        void ft(const std::valarray<double>& u, const std::valarray<double>& v);
        std::pair<double,double> center_mass() const;
        void shift_xy(std::pair<double, double> xy);
        void recenter();
        std::valarray<double> get_mu_real() const { return mu_real; }
        std::valarray<double> get_mu_imag() const { return mu_imag; }
        void print(std::ostream& out) const;
        std::string description() const;
        void from_prior(DNest4::RNG& rng);
        double perturb(DNest4::RNG& rng);

    private:
        std::vector<CGComponent> components_;
        // SkyModel prediction
        std::valarray<double> mu_real;
        std::valarray<double> mu_imag;

};


#endif //BSC_MODEL_H

#ifndef BSC_MODEL_H
#define BSC_MODEL_H

#include <vector>
#include "Component.h"
#include "RNG.h"


class SkyModel {
    public:
        SkyModel();
        SkyModel(const SkyModel& other);
        ~SkyModel();
        SkyModel&operator=(const SkyModel& other);

        void add_component(Component* component);
        void ft_update(const std::valarray<double>& u, const std::valarray<double>& v);
        void ft(const std::valarray<double>& u, const std::valarray<double>& v);
        std::pair<double,double> center_mass() const;
        std::pair<double,double> center_mass2() const;
        void shift_xy(std::pair<double, double> xy);
        void recenter();
        std::valarray<double> get_mu_real() const { return mu_real; }
        std::valarray<double> get_mu_imag() const { return mu_imag; }
        size_t size() const;
        void print(std::ostream& out) const;
        std::string description() const;
        void from_prior(DNest4::RNG& rng);
        // MH proposal for SkyModel. Returns logH
        double perturb(DNest4::RNG& rng);
        int get_n_components();
        void find_brightest();
        void phase_shift_mu(std::pair<double, double> shift);

    private:
        std::vector<Component*> components_;
        std::vector<int> components_sizes_;
        // SkyModel prediction
        std::valarray<double> mu_real;
        std::valarray<double> mu_imag;
        int idx_brightest{-1};
        int idx_brightest_old{-1};
        std::pair<double,double> last_shift{};

};


#endif //BSC_MODEL_H

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
        void ft(const std::valarray<double>& u, const std::valarray<double>& v);
        void ft_from_all(const std::valarray<double>& u, const std::valarray<double>& v);
        std::pair<double,double> center_brightest() const;
        void shift_xy(std::pair<double, double> xy);
        void recenter_brightest();
        void update_brightest();
        void set_param_vec(std::valarray<double> param);
        std::valarray<double> get_mu_real() const { return mu_real; }
        std::valarray<double> get_mu_imag() const { return mu_imag; }
        size_t size() const;
        void print(std::ostream& out) const;
        std::string description() const;
        void from_prior(DNest4::RNG& rng);
        // MH proposal for SkyModel. Returns logH
        double perturb(DNest4::RNG& rng);
        int get_n_components();
        std::vector<int> get_components_sizes();
        void set_x(double x);
        // This is true if position of the brightest component is changed:
        // 1. The position of brightest component changes (x or y)
        // 2. The flux of brightest component decreases and other component becomes the brightest
        // 3. The flux of other component increases and it becomes the brightest
        bool position_of_brightest_is_updated;

    private:
        std::vector<Component*> components_;
        std::vector<int> components_sizes_;
        // SkyModel prediction
        std::valarray<double> mu_real;
        std::valarray<double> mu_imag;
        // This phase shift all components predictions after re-centering
        std::valarray<double> cos_theta;
        std::valarray<double> sin_theta;

};


#endif //BSC_MODEL_H

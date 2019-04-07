#ifndef BSC_MODEL_H
#define BSC_MODEL_H

#include <vector>
#include "Component.h"


class SkyModel {
    public:
        SkyModel();
        void add_component(Component* component);
        void ft(const std::valarray<double>& u, const std::valarray<double>& v);
        void set_param_vec(std::valarray<double> param);
        std::valarray<double> get_mu_real() const { return mu_real; }
        std::valarray<double> get_mu_imag() const { return mu_imag; }
        size_t size() const;
        void print(std::ostream& out) const;

    private:
        std::vector<Component*> components_;
        // SkyModel prediction
        std::valarray<double> mu_real;
        std::valarray<double> mu_imag;

};


#endif //BSC_MODEL_H

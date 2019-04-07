#include <SkyModel.h>
#include <iostream>


SkyModel::SkyModel() = default;

void SkyModel::add_component(Component *component) {
    components_.push_back(component);
}

void SkyModel::ft(const std::valarray<double>& u, const std::valarray<double>& v) {
    std::valarray<double> real (0.0, u.size());
    std::valarray<double> imag (0.0, u.size());
    for (auto comp : components_) {
        comp->ft(u, v);
        real = real + comp->get_mu_real();
        imag = imag + comp->get_mu_imag();
    }
    mu_real = real;
    mu_imag = imag;

}

void SkyModel::set_param_vec(std::valarray<double> param) {
    size_t size_used = 0;

    for (auto comp : components_) {
        comp->set_param_vec(param[std::slice(size_used, comp->size(), 1)]);
        size_used += comp->size();
    }

}

size_t SkyModel::size() const {
    size_t cursize = 0;
    if (components_.empty()) {
        return cursize;
    }
    else {
        for (auto comp : components_) {
            cursize += comp->size();
        }
    }
    return cursize;
}

void SkyModel::print(std::ostream &out) const
{
    for (auto comp: components_) {
        comp->print(out);
    }
}

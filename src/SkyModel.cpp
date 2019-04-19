#include <SkyModel.h>
#include <iostream>


SkyModel::SkyModel() = default;


SkyModel::SkyModel(const SkyModel &other) {
    for (auto other_comp : other.components_) {
        Component* comp = other_comp->clone();
        components_.push_back(comp);
    }
}


void SkyModel::add_component(Component *component) {
    components_.push_back(component);
    components_sizes_.push_back(component->size());
}


void SkyModel::ft(const std::valarray<double>& u, const std::valarray<double>& v) {
    std::valarray<double> real (0.0, u.size());
    std::valarray<double> imag (0.0, u.size());
    for (auto comp : components_) {
        comp->ft(u, v);
        real = real + comp->get_mu_real();
        imag = imag + comp->get_mu_imag();
        std::cout << "Component real[0] = " << real[0] << std::endl;
    }
    mu_real = real;
    mu_imag = imag;
    std::cout << "SkuModel mu_real[0] = " << mu_real[0] << std::endl;

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


void SkyModel::from_prior(DNest4::RNG &rng) {
    std::cout << "Generating from prior SkyModel" << std::endl;
    for (auto comp: components_) {
        comp->from_prior(rng);
    }
}


double SkyModel::perturb(DNest4::RNG &rng) {
    int which = rng.rand_int(get_n_components());
    return components_[which]->perturb(rng);
}


int SkyModel::get_n_components() {
    return components_.size();
}


std::vector<int> SkyModel::get_components_sizes() {
    return components_sizes_;
}


void SkyModel::set_x(double x) {
    components_[0]->set_x(x);
}



#include <SkyModel.h>
#include <iostream>


SkyModel::SkyModel() = default;


SkyModel::SkyModel(const SkyModel &other) {
    for (auto other_comp : other.components_) {
        Component* comp = other_comp->clone();
        components_.push_back(comp);
    }
    components_sizes_ = other.components_sizes_;
    mu_real = other.mu_real;
    mu_imag = other.mu_imag;
}


SkyModel& SkyModel::operator=(const SkyModel& other) {
    if (&other == this) {
        return *this;
    }

    for (auto comp : components_) {
        delete comp;
    }
    components_.clear();

    for (auto other_comp : other.components_) {
        Component* comp = other_comp->clone();
        components_.push_back(comp);
    }
    components_sizes_ = other.components_sizes_;
    mu_real = other.mu_real;
    mu_imag = other.mu_imag;
    return *this;
}



SkyModel::~SkyModel() {
    for (auto comp : components_) {
        delete comp;
    }
    components_.clear();
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
    }
    mu_real = real;
    mu_imag = imag;

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
    for (auto comp: components_) {
        comp->from_prior(rng);
    }
    // Move center of mass to the phase center
    auto xy = center_mass();
    shift_xy(xy);
}


double SkyModel::perturb(DNest4::RNG &rng) {
    int which = rng.rand_int(get_n_components());
    return components_[which]->perturb(rng);
}


int SkyModel::get_n_components() {
    return components_.size();
}



std::string SkyModel::description() const {
    std::string descr;
    for (auto comp: components_) {
        descr += comp->description();
        descr += " ";
    }
    descr.pop_back();
    return descr;
}


std::pair<double, double> SkyModel::center_mass() const {
    double x_b = 0;
    double y_b = 0;
    double flux_b = 0;
    for (auto comp: components_) {
        if(comp->get_flux() > flux_b) {
            x_b = comp->get_x();
            y_b = comp->get_y();
            flux_b = comp->get_flux();
        }
    }
    return std::make_pair<double, double>(reinterpret_cast<double &&>(x_b), reinterpret_cast<double &&>(y_b));
}


void SkyModel::shift_xy(std::pair<double, double> xy) {
    for (auto comp: components_) {
        comp->shift_xy(xy);
    }
}


void SkyModel::recenter() {
    auto xy = center_mass();
    shift_xy(xy);
}



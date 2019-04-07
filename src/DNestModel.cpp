#include "DNestModel.h"


DNestModel::DNestModel() {

}


void DNestModel::from_prior(DNest4::RNG &rng) {

}


double DNestModel::perturb(DNest4::RNG &rng) {
    return 0;
}


double DNestModel::log_likelihood() const {
    return 0;
}


void DNestModel::print(std::ostream &out) const {

}


std::string DNestModel::description() const {
    return std::__cxx11::string();
}

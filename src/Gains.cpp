#include <iostream>
#include "Gains.h"
#include "Gain.h"



//Gains::Gains() = default;


Gains::Gains(Data data) {
    // This map from antenna numbers to antenna indexes.
    auto antennas_map_inv = data.get_antennas_map_inv();
    for (int i=0; i<data.n_antennas(); i++) {
        if (antennas_map_inv[i] == data.constant_gain_antenna) {
            gains.emplace_back(new ConstantGain(data.get_times_amp()[antennas_map_inv[i]], data.get_times_phase()[antennas_map_inv[i]]));

        } else {
            gains.emplace_back(new Gain(data.get_times_amp()[antennas_map_inv[i]], data.get_times_phase()[antennas_map_inv[i]]));
        }
    }
}


Gains::Gains(Gains &other) {
    gains = std::vector<BasicGain*>();
    for (auto gain : other.gains) {
        BasicGain* b = gain->clone();
        gains.emplace_back(b);
    }
}


Gains& Gains::operator=(const Gains& other) {
    if (&other == this)
        return *this;

    for (auto gain : gains) {
        delete gain;
    }

    gains = std::vector<BasicGain*>();
    for (auto gain : other.gains) {
        BasicGain* b = gain->clone();
        gains.emplace_back(b);
    }
    return *this;
}


Gains::~Gains() {
    for (auto gain : gains) {
        delete gain;
    }
}


int Gains::size() const {
    return gains.size();
}


BasicGain* Gains::operator[](int i) {
    return gains[i];
}


void Gains::from_prior(DNest4::RNG &rng) {
    for (auto gain: gains) {
        gain->from_prior(rng);
    }
}


double Gains::perturb(DNest4::RNG& rng) {
    double logH = 0;
    // Leaving first gains as it is (amplitudes = 1s, phases = 0s)
    int which = rng.rand_int(gains.size()-1);
    //std::cout << "Perturbing gain # " << which+1 << std::endl;
    logH += gains[which+1]->perturb(rng);
    return logH;
}


std::string Gains::description() const {
    std::string descr;
    for (auto gain : gains) {
        descr += gain->description();
        descr += " ";
    }
    descr.pop_back();
    return descr;
}


void Gains::print(std::ostream &out) const {
    for (auto gain : gains) {
        gain->print(out);
    }
}



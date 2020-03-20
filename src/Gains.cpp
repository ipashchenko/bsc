#include <iostream>
#include "Gains.h"
#include "Gain.h"
#include "MyExceptions.h"



Gains::Gains(Data data, int refant_) {
    // ant_i of the reference antenna
    refant = refant_;

    // Map between ant_i and its position in antennas vector.
    auto antennas_map = data.get_antennas_map();
    if(!antennas_map.count(refant)) {
        throw BadReferenceAntenna();
    }

    for (int i=0; i<data.n_antennas(); i++) {
        if (i != antennas_map[refant]) {
            antennas_changing_gain.push_back(i);
        }
    }

    // Map between antenna position in antennas vector and ant_i.
    auto antennas_map_inv = data.get_antennas_map_inv();
    for (int i=0; i<data.n_antennas(); i++) {
        gains.emplace_back(new Gain(antennas_map_inv[i], data.get_times_amp()[antennas_map_inv[i]], data.get_times_phase()[antennas_map_inv[i]]));
    }
}


Gains::Gains(Gains &other) {
    gains = std::vector<Gain*>();
    antennas_changing_gain = other.antennas_changing_gain;
    refant = other.refant;
    for (auto gain : other.gains) {
        gains.emplace_back(new Gain(*gain));
    }
}


Gains& Gains::operator=(const Gains& other) {
    if (&other == this)
        return *this;

    for (auto gain : gains) {
        delete gain;
    }

    gains = std::vector<Gain*>();
    for (auto gain : other.gains) {
        gains.emplace_back(new Gain(*gain));
    }
    antennas_changing_gain = other.antennas_changing_gain;
    refant = other.refant;
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


Gain* Gains::operator[](int i) {
    return gains[i];
}


void Gains::from_prior_amp_mean(DNest4::RNG &rng) {
    for (auto i : antennas_changing_gain) {
        gains[i]->from_prior_amp_mean(rng);
    }
}


void Gains::from_prior_phase_mean(DNest4::RNG &rng) {
    for (auto i : antennas_changing_gain) {
        gains[i]->from_prior_phase_mean(rng);
    }
}


void Gains::from_prior_amp(DNest4::RNG &rng) {
    for (auto i : antennas_changing_gain) {
        gains[i]->from_prior_amp(rng);
    }
}


void Gains::from_prior_phase(DNest4::RNG &rng) {
    for (auto i : antennas_changing_gain) {
        gains[i]->from_prior_phase(rng);
    }
}


void Gains::sum() {
    for (auto i : antennas_changing_gain) {
        gains[i]->sum();
    }
}


double Gains::perturb(DNest4::RNG& rng) {
    double logH = 0;
    int which = rng.rand_int(antennas_changing_gain.size());
    logH += gains[antennas_changing_gain[which]]->perturb(rng);
    gains[antennas_changing_gain[which]]->sum();
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



#include <utility>
#include <random>
#include "Gain.h"
#include <Eigen/Core>
#include <Eigen/Eigen>
#include <Gain.h>
#include <valarray>
#include <set>
#include <iostream>
#include <Distributions/Gaussian.h>
#include <Distributions/Uniform.h>
#include "RNG.h"
#include "Utils.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;


Gain::Gain(int ant_, const std::set<double>& new_times_amp, const std::set<double>& new_times_phase) :
    ant(ant_),
    phase_mean(0.0),
    amp_mean(1.0)

    {
    // Get unique times of gain amplitudes
    std::vector<double> times_amp_vec;
    times_amp_vec.assign(new_times_amp.begin(), new_times_amp.end());
    // Get unique times of gain phases
    std::vector<double> times_phase_vec;
    times_phase_vec.assign(new_times_phase.begin(), new_times_phase.end());
    // Convert to Eigen VectorXd
    times_amp = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(times_amp_vec.data(), times_amp_vec.size());
    times_phase = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(times_phase_vec.data(), times_phase_vec.size());

    amplitudes = std::valarray<double>(0.0, times_amp.size());
    phases = std::valarray<double>(0.0, times_amp.size());
    amplitudes_full = std::valarray<double>(0.0, times_amp.size());
    phases_full = std::valarray<double>(0.0, times_amp.size());
    }


Gain::Gain(int ant_, const std::set<double>& new_times) :
    ant(ant_),
    phase_mean(0.0),
    amp_mean(1.0)

    {
    // Get unique times of gain amplitude & phase
    std::vector<double> times_vec;
    times_vec.assign(new_times.begin(), new_times.end());
    // Convert to Eigen VectorXd
    times_amp = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(times_vec.data(), times_vec.size());
    times_phase = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(times_vec.data(), times_vec.size());

    amplitudes = std::valarray<double>(0.0, times_amp.size());
    phases = std::valarray<double>(0.0, times_amp.size());
    amplitudes_full = std::valarray<double>(0.0, times_amp.size());
    phases_full = std::valarray<double>(0.0, times_amp.size());
    }


void Gain::print_times(std::ostream &out) const
{
    out << "times amp = " << std::endl;
    for (int i=0; i<times_amp.size(); i++) {
        out << times_amp[i] << ", ";
    }
    out << std::endl;

    out << "times phase = " << std::endl;
    for (int i=0; i<times_phase.size(); i++) {
        out << times_phase[i] << ", ";
    }
    out << std::endl;

}


void Gain::from_prior_amp(DNest4::RNG &rng) {
    amplitudes = make_normal_random(size_amp(), 0.0, 0.05, rng);
}


void Gain::from_prior_phase(DNest4::RNG &rng) {
    phases = make_normal_random(size_phase(), 0.0, 0.1, rng);
}


void Gain::from_prior_amp_mean(DNest4::RNG& rng) {
    amp_mean = 1.0 + 0.1*rng.rand();
}


void Gain::from_prior_phase_mean(DNest4::RNG& rng) {
    const DNest4::Uniform unif(-1.0*M_PI, 1.0*M_PI);
    phase_mean = unif.generate(rng);
}


int Gain::size_amp() {
    return times_amp.size();
}


int Gain::size_phase() {
    return times_phase.size();
}


void Gain::sum() {
    amplitudes_full = amplitudes + amp_mean;
    phases_full = phases + phase_mean;
}

const std::valarray<double>& Gain::get_amplitudes() const {
    return amplitudes_full;
}


const std::valarray<double>& Gain::get_phases() const {
    return phases_full;
}

void Gain::print_amplitudes(std::ostream &out) const {
    out << "amplitudes = " << std::endl;
    for (int i=0; i < times_amp.size(); i++) {
        out << amplitudes[i] << ", ";
    }
    out << std::endl;
}


void Gain::print_phases(std::ostream &out) const {
    out << "phases = " << std::endl;
    for (int i=0; i < times_phase.size(); i++) {
        out << phases[i] << ", ";
    }
    out << std::endl;
}


double Gain::perturb(DNest4::RNG &rng) {
    double logH = 0.;

    int which = rng.rand_int(2);
    // Amplitude GP
    if(which == 0) {

        // More often perturb phase latent variables
        if(rng.rand() <= 0.9) {
            // Choose what value to perturb
            int which = rng.rand_int(amplitudes.size());

            logH -= -0.5*pow(amplitudes[which]/0.05, 2.0);
            amplitudes[which] += 0.05*rng.randn();
            logH += -0.5*pow(amplitudes[which]/0.05, 2.0);

            // Pre-reject
            if(rng.rand() >= exp(logH)) {
                return -1E300;
            } else
                logH = 0.0;
        } else {
            logH -= -0.5*pow((amp_mean-1.0)/0.1, 2.0);
            amp_mean += 0.1*rng.randh();
            logH += -0.5*pow((amp_mean-1.0)/0.1, 2.0);

            // Pre-reject
            if (rng.rand() >= exp(logH)) {
                return -1E300;
            } else
                logH = 0.0;
        }

    }
    // Phase GP
    else {

        // More often perturb phase latent variables
        if(rng.rand() <= 0.9) {

            // Choose what value to perturb
            int which = rng.rand_int(phases.size());

            logH -= -0.5*pow(phases[which]/0.1, 2.0);
            phases[which] += 0.1*rng.randn();
            logH += -0.5*pow(phases[which]/0.1, 2.0);

            // Pre-reject
            if (rng.rand() >= exp(logH)) {
                return -1E300;
            } else
                logH = 0.0;

        } else {
            phase_mean += 2.0*M_PI*rng.randh();
            DNest4::wrap(phase_mean, -1.0*M_PI, M_PI);

            // Pre-reject
            if (rng.rand() >= exp(logH)) {
                return -1E300;
            } else
                logH = 0.0;
        }

    }

    return logH;
}


std::string Gain::description() const {

    std::string suffix;
    suffix += "_ant" + std::to_string(ant);

    std::string descr;
    descr += "amp_mean" + suffix + " ";
    for (int i = 0; i < times_amp.size(); i++) {
        descr += ("amp" + std::to_string(i) + suffix + " ");
    }
    descr += "phase_mean " + suffix + " ";
    for (int i = 0; i < times_phase.size(); i++) {
        descr += ("phase" + std::to_string(i) + suffix + " ");
    }
    descr.pop_back();
    return descr;
}


void Gain::print(std::ostream &out) const {
    out << amp_mean << '\t';
    for (int i = 0; i < times_amp.size(); i++) {
        out << amplitudes[i] << '\t';
    }
    out << phase_mean << '\t';
    for (int i = 0; i < times_phase.size(); i++) {
        out << phases[i] << '\t';
    }
}


std::valarray<double> make_normal_random(int number, double mu, double scale, DNest4::RNG& rng)
{
    const DNest4::Gaussian gauss(mu, scale);
    std::valarray<double> randNums(number);
    for (int i=0; i<number; i++) {
        randNums[i] = gauss.generate(rng);
    }
    return randNums;
}
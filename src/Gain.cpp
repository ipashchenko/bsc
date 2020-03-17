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


Gain::Gain(int ant_, int IF_, const std::set<double>& new_times_amp, const std::set<double>& new_times_phase) :
    ant(ant_),
    IF(IF_),
    logamp_amp(-3.0),
    logamp_phase(-2.0),
    logscale_amp(5.0),
    logscale_phase(4.0),
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
    v_amp = std::valarray<double>(0.0, times_amp.size());
    v_phase = std::valarray<double>(0.0, times_phase.size());
    C_amp = MatrixXd(times_amp.size(), times_amp.size());
    C_phase = MatrixXd(times_phase.size(), times_phase.size());
    L_amp = MatrixXd(times_amp.size(), times_amp.size());
    L_phase = MatrixXd(times_phase.size(), times_phase.size());

    calculate_C_amp();
    calculate_C_phase();
    calculate_L_amp();
    calculate_L_phase();
    }


Gain::Gain(int ant_, int IF_, const std::set<double>& new_times) :
    ant(ant_),
    IF(IF_),
    logamp_amp(-3.0),
    logamp_phase(-2.0),
    logscale_amp(5.0),
    logscale_phase(5.0),
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
    v_amp = std::valarray<double>(0.0, times_amp.size());
    v_phase = std::valarray<double>(0.0, times_phase.size());
    C_amp = MatrixXd(times_amp.size(), times_amp.size());
    C_phase = MatrixXd(times_phase.size(), times_phase.size());
    L_amp = MatrixXd(times_amp.size(), times_amp.size());
    L_phase = MatrixXd(times_phase.size(), times_phase.size());

    calculate_C_amp();
    calculate_C_phase();
    calculate_L_amp();
    calculate_L_phase();
    }


void Gain::from_prior_v_amp(DNest4::RNG &rng) {
    v_amp = make_normal_random(size_amp(), rng);
}


void Gain::from_prior_v_phase(DNest4::RNG &rng) {
    v_phase = make_normal_random(size_phase(), rng);
}


void Gain::from_prior_amp_mean(DNest4::RNG& rng) {
    amp_mean = 1.0 + 0.1*rng.rand();
}


void Gain::from_prior_phase_mean(DNest4::RNG& rng) {
    const DNest4::Uniform unif(-1.0*M_PI, 1.0*M_PI);
    phase_mean = unif.generate(rng);
}


void Gain::from_prior_hp_amp(DNest4::RNG& rng) {
    //logamp_amp = -3.0 + 1.0*rng.randn();
    //logscale_amp = 7.0 + 1.0*rng.randn();
    logamp_amp = -3.0;
    logscale_amp = 5.0;
}


void Gain::from_prior_hp_phase(DNest4::RNG& rng) {
    //logamp_phase = -3.0 + 1.0*rng.randn();
    //logscale_phase = 5.0 + 1.0*rng.randn();
    logamp_phase = -2.0;
    logscale_phase = 4.0;
}


void Gain::calculate_C_amp() {
    Eigen::MatrixXd sqdist = - 2*times_amp*times_amp.transpose();
    sqdist.rowwise() += times_amp.array().square().transpose().matrix();
    sqdist.colwise() += times_amp.array().square().matrix();
    sqdist *= (-0.5/(exp(2.0*logscale_amp)));
    C_amp = exp(logamp_amp) * sqdist.array().exp();
}


void Gain::calculate_C_phase() {
    Eigen::MatrixXd sqdist = - 2*times_phase*times_phase.transpose();
    sqdist.rowwise() += times_phase.array().square().transpose().matrix();
    sqdist.colwise() += times_phase.array().square().matrix();
    sqdist *= (-0.5/(exp(2.0*logscale_phase)));
    C_phase = exp(logamp_phase) * sqdist.array().exp();
}


void Gain::calculate_L_amp() {
    // perform the Cholesky decomposition of covariance matrix
    Eigen::LLT<Eigen::MatrixXd> cholesky = C_amp.llt();
    // get the lower triangular matrix L
    L_amp = cholesky.matrixL();
}


void Gain::calculate_L_phase() {
    // perform the Cholesky decomposition of covariance matrix
    Eigen::LLT<Eigen::MatrixXd> cholesky = C_phase.llt();
    // get the lower triangular matrix L
    L_phase = cholesky.matrixL();
}


int Gain::size_amp() {
    return times_amp.size();
}


int Gain::size_phase() {
    return times_phase.size();
}


std::valarray<double>& Gain::get_amplitudes() {
    return amplitudes;
}


std::valarray<double>& Gain::get_phases() {
    return phases;
}


void Gain::calculate_amplitudes() {
    // Convert std::valarray ``v_amp`` to Eigen::VectorXd
    VectorXd v_amp_vec = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(&v_amp[0], v_amp.size());
    VectorXd amp = L_amp*v_amp_vec;
    // Convert Eigen::VectorXd to std::valarray
    amplitudes = amp_mean + std::valarray<double>(amp.data(), amp.size());
}


void Gain::calculate_phases() {
    // Convert std::valarray ``v_amp`` to Eigen::VectorXd
    VectorXd v_phase_vec = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(&v_phase[0], v_phase.size());
    VectorXd phase = L_phase*v_phase_vec;
    // Convert Eigen::VectorXd to std::valarray
    phases = phase_mean + std::valarray<double>(phase.data(), phase.size());
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

            logH -= -0.5*pow(v_amp[which], 2.0);
            v_amp[which] += rng.randn();
            logH += -0.5*pow(v_amp[which], 2.0);

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

        calculate_amplitudes();
    }
    // Phase GP
    else {

        // More often perturb phase latent variables
        if(rng.rand() <= 0.9) {

            // Choose what value to perturb
            int which = rng.rand_int(phases.size());

            logH -= -0.5*pow(v_phase[which], 2.0);
            v_phase[which] += rng.randn();
            logH += -0.5*pow(v_phase[which], 2.0);

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

        calculate_phases();

    }

    return logH;
}


std::string Gain::description() const {
    std::string suffix;
    suffix += "_ant" + std::to_string(ant) + "_" + "IF" + std::to_string(IF);
    std::string descr;
    descr += "amp_mean" + suffix + " ";
    for (int i = 0; i < times_amp.size(); i++) {
        descr += ("amp" + std::to_string(i) + suffix + " ");
    }
    descr += "phase_mean ";
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


std::valarray<double> make_normal_random(int number)
{
    std::random_device rd;
    std::normal_distribution<double> normalDistr(0, 1);
    std::minstd_rand generator(rd());
    std::valarray<double> randNums(number);

    for (int i=0; i<number; i++) {
        randNums[i] = normalDistr(generator);
    }
    return randNums;
}


std::valarray<double> make_normal_random(int number, DNest4::RNG& rng)
{
    const DNest4::Gaussian gauss(0.0, 1.0);
    std::valarray<double> randNums(number);
    for (int i=0; i<number; i++) {
        randNums[i] = gauss.generate(rng);
    }
    return randNums;
}
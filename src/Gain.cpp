#include <utility>
#include <random>
#include "Gain.h"
#include <Eigen/Core>
#include <Eigen/Eigen>
#include <Gain.h>
#include <valarray>
#include <set>
#include <iostream>
#include "RNG.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;


Gain::Gain(std::vector<double> new_times_amp, std::vector<double> new_times_phase) :
    logamp_amp(-3.0),
    logamp_phase(-3.0),
    logscale_amp(7.0),
    logscale_phase(5.0)

    {
    // Get unique times of gain amplitudes
    std::vector<double> times_amp_vec;
    std::set<double> s_amp( new_times_amp.begin(), new_times_amp.end() );
    times_amp_vec.assign( s_amp.begin(), s_amp.end() );
    // Get unique times of gain phases
    std::vector<double> times_phase_vec;
    std::set<double> s_phase( new_times_phase.begin(), new_times_phase.end() );
    times_phase_vec.assign( s_phase.begin(), s_phase.end() );
    // Convert to Eigen VectorXd
    times_amp = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(times_amp_vec.data(), times_amp_vec.size());
    times_phase = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(times_phase_vec.data(), times_phase_vec.size());

    amplitudes = std::valarray<double>(1.0, times_amp.size());
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


Gain::Gain(std::vector<double> new_times) :
    logamp_amp(-3.0),
    logamp_phase(-3.0),
    logscale_amp(7.0),
    logscale_phase(5.0)

    {
    // Get unique times of gain amplitude & phase
    std::vector<double> times_vec;
    std::set<double> s( new_times.begin(), new_times.end() );
    times_vec.assign( s.begin(), s.end() );
    // Convert to Eigen VectorXd
    times_amp = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(times_vec.data(), times_vec.size());
    times_phase = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(times_vec.data(), times_vec.size());

    amplitudes = std::valarray<double>(1.0, times_amp.size());
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


//Gain::Gain(const Gain &other) {
//    logamp_amp = other.logamp_amp;
//    logamp_phase = other.logamp_phase;
//    logscale_amp = other.logscale_amp;
//    logscale_phase = other.logscale_phase;
//    times_amp = other.times_amp;
//    times_phase = other.times_phase;
//    amplitudes =other.amplitudes;
//    phases = other.phases;
//    v_amp = other.v_amp;
//    v_phase = other.v_phase;
//    C_amp = other.C_amp;
//    C_phase = other.C_phase;
//    L_amp = other.L_amp;
//    L_phase = other.L_phase;
//}


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


void Gain::print_hp(std::ostream &out) const
{
    out << "For amplitude: logamplitude = " << logamp_amp << ", " << "logscale = " << logscale_amp << std::endl;
    out << "For phase: logamplitude = " << logamp_phase << ", " << "logscale = " << logscale_phase << std::endl;

}


void Gain::print_v(std::ostream &out) const
{
    out << "v amp = " << std::endl;
    for (int i=0; i<v_amp.size(); i++) {
        out << v_amp[i] << ", ";
    }
    out << std::endl;

    out << "v_phase = " << std::endl;
    for (int i=0; i<v_phase.size(); i++) {
        out << v_phase[i] << ", ";
    }
    out << std::endl;
}


void Gain::print_C(std::ostream &out) const
{
    out << "C amp = " << std::endl;
    out << C_amp << std::endl;
    out << "C phase = " << std::endl;
    out << C_phase << std::endl;
}


void Gain::print_L(std::ostream &out) const
{
    out << "L amp = " << std::endl;
    out << L_amp << std::endl;
    out << "L phase= " << std::endl;
    out << L_phase << std::endl;
}


void Gain::set_hp_amp(std::valarray<double> params) {
    logamp_amp = params[0];
    logscale_amp = params[1];
}


void Gain::set_hp_phase(std::valarray<double> params) {
    logamp_phase = params[0];
    logscale_phase = params[1];
}


//void Gain::set_v_amp(std::valarray<double> params) {
//    v_amp = std::move(params);
//}
//
//
//void Gain::set_v_phase(std::valarray<double> params) {
//    v_phase = std::move(params);
//}


void Gain::from_prior_v_amp() {
    //v_amp = make_normal_random(size_amp());
    v_amp = std::valarray<double> (0.0, size_amp());
}


void Gain::from_prior_v_phase() {
    //v_phase = make_normal_random(size_phase());
    v_phase = std::valarray<double> (0.0, size_phase());
}


void Gain::from_prior_hp_amp(DNest4::RNG& rng) {
    //std::cout << "Generating from prior Gain HP AMP" << std::endl;
    //logamp_amp = -3.0 + 1.0*rng.randn();
    //logscale_amp = 7.0 + 1.0*rng.randn();
    logamp_amp = -3.0;
    logscale_amp = 7.0;
}


void Gain::from_prior_hp_phase(DNest4::RNG& rng) {
    //logamp_phase = -3.0 + 1.0*rng.randn();
    //logscale_phase = 5.0 + 1.0*rng.randn();
    logamp_phase = -3.0;
    logscale_phase = 5.0;
}


void Gain::calculate_C_amp() {
    //std::cout << "Calculating C of amp GP" << std::endl;
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
    //std::cout << "Calculating L of amp GP" << std::endl;
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
    //std::cout << "Calculating amplitudes" << std::endl;
    // Convert std::valarray ``v_amp`` to Eigen::VectorXd
    VectorXd v_amp_vec = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(&v_amp[0], v_amp.size());


    //// Debug print out v
    //std::valarray<double> v_amp = std::valarray<double>(v_amp_vec.data(), v_amp_vec.size());
    ////std::cout << "V amplitudes = " << std::endl;
    //for (int i=0; i < times_amp.size(); i++) {
    //    std::cout << v_amp[i] << ", ";
    //}
    //std::cout << std::endl;


    VectorXd amp = L_amp*v_amp_vec;
    // Convert Eigen::VectorXd to std::valarray
    amplitudes = 1.0 + std::valarray<double>(amp.data(), amp.size());
    //std::cout << "Calculated amplitudes: " << std::endl;
    //print_amplitudes(std::cout);

}


void Gain::calculate_phases() {
    // Convert std::valarray ``v_amp`` to Eigen::VectorXd
    VectorXd v_phase_vec = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(&v_phase[0], v_phase.size());

    VectorXd phase = L_phase*v_phase_vec;


    //// Debug print out v
    //std::valarray<double> v_phases = std::valarray<double>(v_phase_vec.data(), v_phase_vec.size());
    //std::cout << "V phases = " << std::endl;
    //for (int i=0; i < times_phase.size(); i++) {
    //    std::cout << v_phases[i] << ", ";
    //}
    //std::cout << std::endl;


    // Convert Eigen::VectorXd to std::valarray
    phases = std::valarray<double>(phase.data(), phase.size());
    //std::cout << "Calculated phases: " << std::endl;
    //print_phases(std::cout);

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

        // Choose what value to perturb
        int which = rng.rand_int(amplitudes.size());

        logH -= -0.5*pow(v_amp[which], 2.0);
        v_amp[which] += rng.randn();
        logH += -0.5*pow(v_amp[which], 2.0);

        // Pre-reject
        if(rng.rand() >= exp(logH)) {
            return -1E300;
        }
        else
            logH = 0.0;

        calculate_amplitudes();
    }
    // Phase GP
    else {

        // Choose what value to perturb
        int which = rng.rand_int(phases.size());

        logH -= -0.5*pow(v_phase[which], 2.0);
        v_phase[which] += rng.randn();
        logH += -0.5*pow(v_phase[which], 2.0);

        // Pre-reject
        if(rng.rand() >= exp(logH)) {
            return -1E300;
        }
        else
            logH = 0.0;

        calculate_phases();

    }

    return logH;
}


std::string Gain::description() const {
    std::string descr;
    for (int i = 0; i < times_amp.size(); i++) {
        descr += ("amp" + std::to_string(i) + " ");
    }
    for (int i = 0; i < times_phase.size(); i++) {
        descr += ("phase" + std::to_string(i) + " ");
    }
    descr.pop_back();
    return descr;
}


void Gain::print(std::ostream &out) const {
    for (int i = 0; i < times_amp.size(); i++) {
        out << amplitudes[i] << '\t';
    }
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
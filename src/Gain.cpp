#include <utility>
#include <random>
#include "Gain.h"
#include <Eigen/Core>
#include <Eigen/Eigen>
#include <Gain.h>
#include <valarray>
#include <set>
#include <iostream>

using Eigen::VectorXd;
using Eigen::MatrixXd;


Gain::Gain(std::vector<double> new_times_amp, std::vector<double> new_times_phase) :
    amp_amp(0.1),
    amp_phase(0.1),
    scale_amp(300),
    scale_phase(30)

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
    amp_amp(1),
    amp_phase(1),
    scale_amp(1),
    scale_phase(1)

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
    out << "For amplitude: amplitude = " << amp_amp << ", " << "scale = " << scale_amp << std::endl;
    out << "For phase: amplitude = " << amp_phase << ", " << "scale = " << scale_phase << std::endl;

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
    amp_amp = params[0];
    scale_amp = params[1];
}


void Gain::set_hp_phase(std::valarray<double> params) {
    amp_phase = params[0];
    scale_phase = params[1];
}


void Gain::set_v_amp(std::valarray<double> params) {
    v_amp = std::move(params);
}


void Gain::set_v_phase(std::valarray<double> params) {
    v_phase = std::move(params);
}


void Gain::from_prior_v_amp() {
    v_amp = make_normal_random(size_amp());
}


void Gain::from_prior_v_phase() {
    v_phase = make_normal_random(size_phase());
}


void Gain::calculate_C_amp() {
    Eigen::MatrixXd sqdist = - 2*times_amp*times_amp.transpose();
    sqdist.rowwise() += times_amp.array().square().transpose().matrix();
    sqdist.colwise() += times_amp.array().square().matrix();
    sqdist *= (-0.5/(scale_amp*scale_amp));
    C_amp = amp_amp * sqdist.array().exp();
}


void Gain::calculate_C_phase() {
    Eigen::MatrixXd sqdist = - 2*times_phase*times_phase.transpose();
    sqdist.rowwise() += times_phase.array().square().transpose().matrix();
    sqdist.colwise() += times_phase.array().square().matrix();
    sqdist *= (-0.5/(scale_phase*scale_phase));
    C_phase = amp_phase * sqdist.array().exp();
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
    amplitudes = std::valarray<double>(amp.data(), amp.size());
}


void Gain::calculate_phases() {
    // Convert std::valarray ``v_amp`` to Eigen::VectorXd
    VectorXd v_phase_vec = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(&v_phase[0], v_phase.size());

    VectorXd phase = L_phase*v_phase_vec;
    // Convert Eigen::VectorXd to std::valarray
    phases = std::valarray<double>(phase.data(), phase.size());
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
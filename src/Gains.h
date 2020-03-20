#ifndef BSC__GAINS_H_
#define BSC__GAINS_H_

#include "Gain.h"
#include "Data.h"
#include "RNG.h"


class Gains {
    public:
        Gains(Data data, int refant);
        Gains(Gains& other);
        Gains& operator=(const Gains& other);
        ~Gains();
        // Generate mean amplitudes and phases
        void from_prior_amp_mean(DNest4::RNG &rng);
        void from_prior_phase_mean(DNest4::RNG &rng);
        void from_prior_amp(DNest4::RNG &rng);
        void from_prior_phase(DNest4::RNG &rng);
        void sum();
        // MH proposal that returns logH
        double perturb(DNest4::RNG& rng);
        int size() const;
        // Getting i-th gain
        Gain* operator[](int i);

        std::string description() const;
        void print(std::ostream &out) const;
    private:
        std::vector<Gain*> gains;
        int refant{};
        std::vector<int> antennas_changing_gain;
};

#endif //BSC__GAINS_H_

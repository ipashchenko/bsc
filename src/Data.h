#ifndef BSC_DATA_H
#define BSC_DATA_H

#include <valarray>
#include <vector>
#include <set>
#include <unordered_map>


class Data
{
    private:
        // Static "global" instance
        static Data instance;
        std::vector<int> antennas;
        // Map between ant_i and its position in antennas vector
        std::unordered_map<int, int> antennas_map;
        // Map between antenna position in antennas vector and ant_i
        std::unordered_map<int, int> antennas_map_inv;
        // Map between ant_i and amplitude times of ant_i
        std::unordered_map<int, std::set<double>> antennas_amp_times_map;
        // Map between ant_i and phases times of ant_i
        std::unordered_map<int, std::set<double>> antennas_phase_times_map;

        // Antenna numbers
        std::vector<int> ant_i;
        std::vector<int> ant_j;
        // uv-coordinates
        std::valarray<double> u;
        std::valarray<double> v;
        // visibilities
        std::valarray<double> vis_real;
        std::valarray<double> vis_imag;
        // error of real/imag part
        std::valarray<double> sigma;
        // original times of visibility measurements
        std::vector<double> times;
        // times of the measurements of gain's amplitude
        std::vector<double> times_amp;
        // Indexes of the antennas ``ant_1`` and ``ant_2`` in the ``times_amp``
        std::vector<int> idx_amp_ant_i;
        std::vector<int> idx_amp_ant_j;
        // times of the measurements of gain's phase
        std::vector<double> times_phase;
        // Indexes of the antennas ``ant_1`` and ``ant_2`` in the ``times_phase``
        std::vector<int> idx_phase_ant_i;
        std::vector<int> idx_phase_ant_j;


    public:
        // Constructor
        Data();
        // Getter for the global instance
        static Data& get_instance()
        { return instance; }
        int n_antennas() const
        { return antennas.size(); }

        // Load data from a file
        void load(const char* filename);

        // Access to the data points
        const std::vector<int>& get_ant_i() const
        { return ant_i; }
        const std::vector<int>& get_ant_j() const
        { return ant_j; }
        const std::valarray<double>& get_u() const
        { return u; }
        const std::valarray<double>& get_v() const
        { return v; }
        const std::valarray<double>& get_vis_real() const
        { return vis_real; }
        const std::valarray<double>& get_vis_imag() const
        { return vis_imag; }
        const std::valarray<double>& get_sigma() const
        { return sigma; }
        const std::vector<double>& get_times() const
        { return times; }

        std::unordered_map<int, std::set<double>>& get_times_amp();

        std::unordered_map<int, std::set<double>>& get_times_phase();

        const std::vector<int>& get_idx_amp_ant_i() const
        { return idx_amp_ant_i; }
        const std::vector<int>& get_idx_amp_ant_j() const
        { return idx_amp_ant_j; }
        const std::vector<int>& get_idx_phase_ant_i() const
        { return idx_phase_ant_i; }
        const std::vector<int>& get_idx_phase_ant_j() const
        { return idx_phase_ant_j; }
        const std::vector<int>& get_antennas() const
        { return antennas; }
        std::unordered_map<int, int>& get_antennas_map()
        { return antennas_map; }
        std::unordered_map<int, int>& get_antennas_map_inv()
        { return antennas_map_inv; }

};


#endif //BSC_DATA_H

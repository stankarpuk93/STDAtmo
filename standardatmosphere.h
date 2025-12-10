#ifndef STANDARDATMOSPHERE_H
#define STANDARDATMOSPHERE_H

#include <array>

namespace standardAtmosphere
{

    constexpr std::array<double, 8> Href = {0, 11000, 20000, 32000, 47000, 51000, 71000, 84852};
    constexpr std::array<double, 8> Kref = {-0.0065, 0, 0.001, 0.0028, 0, -0.0028, -0.002, -0.002};
    constexpr std::array<double, 8> Tiso = {288.15, 216.66, 216.66, 228.65, 270.65, 270.65, 214.65, 187.65};
    constexpr std::array<double, 8> piso = {101325, 22632, 5474, 868.1, 110.9, 66.94, 3.956, 0.373};

    const double a0    = 340.294;
    const double ro0   = 1.225;
    const double p0    = 101325;
    const double g0    = 9.807;
    const double R     = 287.04;
    const double mu0   = 17.16E-6;
    const double S     = 110.6;
    const double T0    = 273.15;
    const double r0    = 6356.766;
    const double gamma = 1.4;
    const double Heps  = 2.0;

    void compute_speed_of_sound_and_viscosity(double &a, double &mu, double T);

    double compute_acceleration_of_gravity(double Hgeo);

    double compute_geopotential_altitude(double Hgeo);

    double compute_geometric_altitude(double Hgp);

    void compute_standard_atmosphere(double H, double dISA, double &g, double &p,
                                     double &T, double &rho, double &mu, double &a);

    double EASfromTAS(double ro, double VTAS);

    double TASfromEAS(double ro, double VEAS);

    double CASfromEAS(double p, double VEAS, double Mach);

    double TASfromMach(double a, double Mach);

    double MachfromTAS(double a, double VTAS);

    double MachfromCAS(double p, double VCAS);

};

#endif // STANDARDATMOSPHERE_H

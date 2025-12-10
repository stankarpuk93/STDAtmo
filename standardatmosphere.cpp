#include "standardatmosphere.h"
#include "cmath"

void standardAtmosphere::compute_speed_of_sound_and_viscosity(double &a, double &mu, double T) {

    // computes speed of sound and dynamic viscosity
    a  = sqrt(gamma * R * T);
    mu = mu0 * pow(T/T0,1.5)*(T0+S)/(T+S);
}

double standardAtmosphere::compute_acceleration_of_gravity(double Hgeo){

    // computes acceleration of gravity wrt altitude
    double g1 = g0 * pow(r0/(r0 + Hgeo / 1000), 2);

    return g1;
}


double standardAtmosphere::compute_geopotential_altitude(double Hgeo){

    // computes a geopotential altitude for a given geometric altitude
    double Hgp = (r0* Hgeo)/(r0 + Hgeo / 1000);

    return Hgp;
}

double standardAtmosphere::compute_geometric_altitude(double Hgp){

    // computes a geometric altitude for a given geopotential altitude
    double Hgeo = (Hgp * r0) / (r0 - Hgp/1000.0);;

    return Hgeo;

}

void standardAtmosphere::compute_standard_atmosphere(double H, double dISA, double &g, double &p,
                                                     double &T, double &rho, double &mu, double &a)
{

    // define all necessary constants
    double rho_iso;
    int Htest_ind = 1;

    // find the index of a required altitude
    while (H > Href[Htest_ind])
    {
        Htest_ind += 1;
    }
    Htest_ind -= 1;

    double K  = Kref[Htest_ind];
    double H1 = Href[Htest_ind];
    double T1 = Tiso[Htest_ind];

    if (K == 0)
    {
        // isothermal region
        double alpha = exp(-g0*(H - H1) / R / T1);
        T   = T1 + dISA;
        p   = piso[Htest_ind] * alpha;
        rho_iso = piso[Htest_ind] / R / T;
        rho = rho_iso * alpha;

    }
    else
    {
        // gradient region
        double alpha = -g0 / R / K;
        double T_noISA = T1 + K * (H - H1);
        T = T_noISA + dISA;
        p = piso[Htest_ind] * pow(T_noISA / T1, alpha);
        rho = p / R / T;
    }

    standardAtmosphere::compute_speed_of_sound_and_viscosity(a, mu, T);
    g = standardAtmosphere::compute_acceleration_of_gravity(H);
}


double standardAtmosphere::EASfromTAS(double ro, double VTAS){

    // converts from TAS to EAS
    double sigma = ro / ro0;
    double VEAS = VTAS * sqrt(sigma);

    return VEAS;
}

double standardAtmosphere::TASfromEAS(double ro, double VEAS){

    // converts from TAS to EAS
    double sigma = ro / ro0;
    double VTAS = VEAS / sqrt(sigma);

    return VTAS;
}

double standardAtmosphere::CASfromEAS(double p, double VEAS, double Mach){

    // converts from EAS to CAS for a given Mach
    double delta = p / p0;
    double qc = p * (pow(1+0.2*Mach*Mach,3.5)-1);
    double VCAS = VEAS * sqrt(1/delta) * sqrt((pow(qc/p0+1,0.2857)-1)/(pow(qc/p+1,0.2857)-1));

    return VCAS;
}

double standardAtmosphere::TASfromMach(double a, double Mach){

    // converts from Mach to TAS
    double VTAS = Mach * a;

    return VTAS;
}

double standardAtmosphere::MachfromTAS(double a, double VTAS){

    // converts from TAS to Mach
    double Mach = VTAS / a;

    return Mach;
}

double standardAtmosphere::MachfromCAS(double p, double CAS){

    // converts from CAS to Mach
    double delta = p / p0;
    double A = pow(1 + 0.5*(gamma-1)*pow(CAS/a0,2),gamma/(gamma-1))-1;
    double Mach = sqrt(2/(gamma-1)*(pow(1/delta*A+1,(gamma-1)/gamma)-1));

    return Mach;
}


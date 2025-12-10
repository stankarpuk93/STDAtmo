#ifndef UNITS_H
#define UNITS_H

class Units
{
public:

    // ========== LENGTH CONVERSIONS (Meters as pivot) ==========
    static double metersToFeet(double meters) { return meters * 3.280839895; }
    static double feetToMeter(double meters) { return meters / 3.280839895; }

    // ========== TEMPERATURE CONVERSIONS (Kelvin as pivot) ==========
    static double kelvinToCelsius(double kelvin) { return kelvin - 273.15; }
    static double kelvinToFahrenheit(double kelvin) { return (kelvin - 273.15) * 9.0/5.0 + 32.0; }
    static double kelvinToRankine(double kelvin) { return kelvin * 9.0/5.0; }

    // ========== PRESSURE CONVERSIONS (Pascals as pivot) ==========
    static double pascalsToPsi(double pascals) { return pascals * 0.00014503773773; }
    static double psiToPascals(double psi) { return psi / 0.00014503773773; }
    static double pascalsToMmHg(double pascals) { return pascals * 0.00750061682704; }
    static double mmHgToPascals(double mmHg) { return mmHg / 0.00750061682704; }

    // ========== DENSITY CONVERSIONS (kg/m³ as pivot) ==========
    static double kgPerM3ToSlugPerFt3(double kgm3) { return kgm3 * 0.00194032033; }
    static double slugPerFt3ToKgPerM3(double slugft3) { return slugft3 / 0.00194032033; }

    // ========== VELOCITY CONVERSIONS (m/s as pivot) ==========
    // Core conversions
    static double metersPerSecondToFeetPerSecond(double ms) { return ms * 3.280839895; }
    static double metersPerSecondToKilometersPerHour(double ms) { return ms * 3.6; }
    static double metersPerSecondToMilesPerHour(double ms) { return ms * 2.2369362921; }
    static double metersPerSecondToKnots(double ms) { return ms * 1.9438444924; }

    static double kilometersPerHourToMetersPerSecond(double kph) { return kph / 3.6; }
    static double feetPerSecondToMetersPerSecond(double fps) { return fps / 3.280839895; }
    static double milesPerHourToMetersPerSecond(double mph) { return mph / 2.2369362921; }
    static double KnotsToMetersPerSecond(double kts) { return kts / 1.9438444924; }


    // ========== VISCOSITY CONVERSIONS (Pa·s as pivot) ==========
    static double pascalSecondToLbfSecondPerFt2(double pas) { return pas * 0.0208854342; }
    static double lbfSecondPerFt2ToPascalSecond(double lbfsft2) { return lbfsft2 / 0.0208854342; }
};

#endif // UNITS_H

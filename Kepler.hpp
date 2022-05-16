//
// Created by artem on 18.04.22.
//
#include <iostream>
#include <cmath>
#include <vector>

class Position {
private:
    double N = 0.0; //Longitude of the ascending node
    double i = 0.0; //Inclination
    double w = 0.0; //Argument of periapsis
    double a = 0.0; //Semimajor axis
    double e = 0.0; //Eccentricity
    double M = 0.0; //Mean anomaly

    double x = 0.0;
    double y = 0.0;
    double z = 0.0;
    double vx = 0.0;
    double vy = 0.0;
    double vz = 0.0;

    //constants
    const double m = 1.32712440019E+20;
    const double AU = 149597870700.0;
    const double tolerance = 0.0000001;
public:

    Position();

    ~Position();

    void set_cartesian(double x, double y, double z, double vx, double vy, double vz);

    void set_kepler(double n, double i, double w, double a, double e, double m);

    double calcE(const double M, const double e);

    void calculating_cartesian();

    void calculating_kepler();

    std::vector<double> get_all_kepler();

    std::vector<double> get_all_cartesian();

};

//
// Created by artem on 18.04.22.
//
#include "Kepler.hpp"


//template<typename T>
class Position {
private:
    double N; //Longitude of the ascending node
    double i; //Inclination
    double w; //Argument of periapsis
    double a; //Semimajor axis
    double e; //Eccentricity
    double M; //Mean anomaly

    double x;
    double y;
    double z;
    double vx;
    double vy;
    double vz;

    double m = 1.32712440019E+20;
public:

    Position(double d) {
        N = rad(ostatok(0.0, 360.0));
        i = rad(0.0);
        w = rad(282.9404 + 4.70935E-5 * d);
        a = 1.000000 * 149597870700;
        e = 0.016709 - 1.151E-9 * d;
        M = rad(ostatok(356.0470 + 0.9856002585 * d, 360.0));
    };

    ~Position() {};

    double rad(double deg) {
        return deg / 180.0 * M_PI;
    }

    double ostatok(const double &delimoe, const double &delitel) {
        double res = delimoe;
        while (res > delitel) {
            res -= delitel;
        }
        return res;
    }

    double calcE(const double M, double e) {
        double E = M;
        while (std::abs(E - M - e * sin(E)) > 0.0001)
            E = E - (E - e * sin(E) - M) / (1 - e * cos(E));
        std::cout << "M = " << M << "   E - e*sinE= " << E - e * sin(E) << "\n";
        return E;
    }

    void Calculating_Cartesian() {

        double E = calcE(M, e);
        double v = 2 * atan2(sqrt(1 + e) * sin(E / 2), sqrt(1 - e) * cos(E / 2));
        double r = a * (1 - e * cos(E));

        double ox = r * cos(v);
        double oy = r * sin(v);
        double oz = 0;
        double ovx = -sqrt(a) / r * sin(E);
        double ovy = sqrt(a) / r * sqrt(1 - e * e) * cos(E);
        double ovz = 0;

        x = ox * (cos(w) * cos(N) - sin(w) * cos(i) * sin(N)) - oy * (sin(w) * cos(N) + cos(w) * cos(i) * sin(N));
        y = ox * (cos(w) * sin(N) + sin(w) * cos(i) * cos(N)) + oy * (cos(w) * cos(i) * cos(N) - sin(w) * sin(N));
        z = ox * (sin(w) * sin(i)) + oy * (cos(w) * sin(i));

        vx = ovx * (cos(w) * cos(N) - sin(w) * cos(i) * sin(N)) - ovy * (sin(w) * cos(w) + cos(w) * cos(i) * sin(N));
        vy = ovx * (cos(w) * sin(w) + sin(w) * cos(i) * cos(N)) + ovy * (cos(w) * cos(i) * cos(w) - sin(w) * sin(N));
        vz = ovx * (sin(w) * sin(i)) + ovy * (cos(w) * sin(i));

    };

    void Calculating_Kepler() {};

    void get_all_Kepler() {
        std::cout << "KEPLER ELEMENTS\n" << "N " << N << "\n" << "i " << (i) << "\n" << "w " << w << "\n" << "a " << a
                  << "\n" << "e " << e
                  << "\n"
                  << "M " << M << "\n";
    };

    void get_all_Cartesian() {
        std::cout << "CARTESIAN COORDINATES\n" << "x " << x << "\n" << "y " << y << "\n" << "z " << z << "\n" << "vx "
                  << vx << "\n" << "vy " << vy << "\n" "vz " << vz << "\n";
    };

};



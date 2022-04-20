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
    double r;
    double v;
    double xh;
    double yh;
    double zh;
    double lonecl;
    double latecl;

public:
    Position() {};

    Position(double d, double ecl) {
        N = rad(ostatok(0.0, 360.0));
        i = rad(0.0);
        w = rad(282.9404 + 4.70935E-5 * d);
        a = 1.000000;// * 149597870700;
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

    void Calculating_Cartesian() {

        double E = M + e * sin(M) * (1.0 + e * cos(M)); //the eccentric anomaly E
        double xv = a * (cos(E) - e);
        double yv = a * (sqrt(1.0 - e * e) * sin(E));

        v = atan2(yv, xv); //the true anomaly v
        r = sqrt(xv * xv + yv * yv);  //the Sun's distance r

        xh = r * (cos(N) * cos(v + w) - sin(N) * sin(v + w) * cos(i));
        yh = r * (sin(N) * cos(v + w) + cos(N) * sin(v + w) * cos(i));
        zh = r * (sin(v + w) * sin(i));

        lonecl = atan2(yh, xh); //the ecliptic longitude
        latecl = atan2(zh, sqrt(xh * xh + yh * yh)); //the ecliptic latitude


    };

    void Calculating_Kepler() {};

    void get_all_Kepler() {
        std::cout << "KEPLER ELEMENTS\n" << "N " << N << "\n" << "i " << i << "\n" << "w " << w << "\n" << "a " << a
                  << "\n" << "e " << e
                  << "\n"
                  << "M " << M << "\n";
    };

    void get_all_Cartesian() {
        std::cout << "CARTESIAN COORDINATES\n" << "r " << r << "\n" << "v " << v << "\n" << "xh " << xh << "\n" << "yh "
                  << yh << "\n" << "zh " << zh << "\n" "lonec1 " << lonecl << "\n" << "latecl " << latecl << "\n";
    };

};



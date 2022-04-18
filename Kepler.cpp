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
    double lonecl;
    double latecl;

public:
    Position() {};

    Position(double d, double ecl) {
        N = rad(0.0);
        i = rad(0.0);
        w = rad(282.9404 + 4.70935E-5 * d);
        a = 1.000000;
        e = 0.016709 - 1.151E-9 * d;
        M = rad(356.0470 + 0.9856002585 * d);
    };

    ~Position() {};

    double rad(double deg) {
        return deg / 180.0 * M_PI;
    }

    void Calculating_Cartesian() {

        double E = M + e * sin(M) * (1.0 + e * cos(M)); //the eccentric anomaly E
        double xv = a * (cos(E) - e);
        double yv = a * (sqrt(1.0 - e * e) * sin(E));

        double v = atan2(yv, xv); //the true anomaly v
        double r = sqrt(xv * xv + yv * yv);  //the Sun's distance r

        double xh = r * (cos(N) * cos(v + w) - sin(N) * sin(v + w) * cos(i));
        double yh = r * (sin(N) * cos(v + w) + cos(N) * sin(v + w) * cos(i));
        double zh = r * (sin(v + w) * sin(i));

        double lonecl = atan2(yh, xh); //the Sun's Right Ascension
        double latecl = atan2(zh, sqrt(xh * xh + yh * yh)); //Declination


    };

    void Calculating_Kepler() {};

    void get_all_Kepler() {
        std::cout << "N " << N << "\n" << "i " << i << "\n" << "w " << w << "\n" << "a " << a << "\n" << "e " << e
                  << "\n"
                  << "M " << M << "\n";
    };

    void get_all_Cartesian() {
        std::cout << "r " << r << "\n" << "v " << v << "\n" << "xh " << xh << "\n" << "yh " << yh << "\n" << "lonec1 "
                  << lonecl << "\n"
                  << "latecl " << latecl << "\n";
    };

};

int main() {
    //DMY on 12th of April 2022 on 9.00
    const double d =
            static_cast<double>(367 * 2022 - 7 * (2022 + (4 + 9) / 12) / 4 + 275 * 4 / 9 + 12 - 730530) +
            9.0 / 24.0;
    const double ecl = 23.4393 - 3.563E-7 * d;

    std::cout << d << "\t" << ecl << "\n";
    Position p = Position(d, ecl);
    p.get_all_Kepler();
    p.Calculating_Cartesian();
    return 0;
}
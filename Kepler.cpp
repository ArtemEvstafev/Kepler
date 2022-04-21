//
// Created by artem on 18.04.22.
//
#include "Kepler.hpp"
#include "overloads.hpp"

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
    Position(double x, double y, double z, double vx, double vy, double vz) :
            x(x), y(y), z(z), vx(vx), vy(vy), vz(vz) {};

    Position(double d) {
        //sun
        /*N = rad(ostatok(0.0, 360.0));
        i = rad(0.0);
        w = rad(282.9404 + 4.70935E-5 * d);
        a = 1.000000 * 149597870700;
        e = 0.016709 - 1.151E-9 * d;
        M = rad(ostatok(356.0470 + 0.9856002585 * d, 360.0));*/
        //mercury
        N = rad(ostatok(48.3313 + 3.24587E-5 * d, 360.0));
        i = rad(7.0047 + 5.00E-8 * d);
        w = rad(29.1241 + 1.01444E-5 * d);
        a = 0.387098 * 149597870700;
        e = 0.205635 + 5.59E-10 * d;
        M = rad(ostatok(168.6562 + 4.0923344368 * d, 360.0));
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

    template<typename T>
    std::vector<T> vectmul(const std::vector<T> &c, const std::vector<T> &b) {
        std::vector<T> result(c.size());
        result[0] = c[1] * b[2] - b[1] * c[2];
        result[1] = c[0] * b[2] - b[0] * c[2];
        result[2] = c[0] * b[1] - b[0] * c[1];
        return result;
    }


    void Calculating_Kepler() {
        std::vector<double> r = {x, y, z};
        std::vector<double> dr = {vx, vy, vz};
        std::vector<double> h = vectmul(r, dr);
        std::cout << "vec h: " << h[0] << " " << h[1] << " " << h[2] << "\n";
        std::vector<double> ve = 1 / m * vectmul(dr, h) - 1 / norm(r) * r;
        std::cout << "vec e: " << ve[0] << " " << ve[1] << " " << ve[2] << "\n";
        std::vector<double> n = {-h[1], h[0], 0.0};
        std::cout << "vec n: " << n[0] << " " << n[1] << " " << n[2] << "\n";
        double v;

        if (r * dr >= 0) {
            if (norm(r) > 0 && norm(dr) > 0)
                v = acos(ve * r / norm(ve) / norm(r));
            else
                v = M_PI / 2.0;
        } else
            v = 2 * M_PI - acos(ve * r / norm(ve) / norm(r));

        if (norm(h) != 0) {
            i = acos(h[2] / norm(h));
        } else
            i = M_PI / 2.0;
        e = norm(ve);
        double E;
        if (e != 0)
            E = atan(tan(v / 2) / sqrt((1 + e) / (1 - e)));
        else
            E = 0;

        if (n[1] >= 0) {
            if (norm(n) > 0)
                N = acos(n[0] / norm(n));
            else
                N = M_PI / 2.0;
        } else
            N = 2 * M_PI - acos(n[0] / norm(n));

        if (ve[2] >= 0) {
            if (norm(n) > 0 && norm(ve) > 0)
                w = acos(ve * n / norm(ve) / norm(n));
            else
                w = M_PI / 2.0;
        } else
            w = 2 * M_PI - acos(ve * n / norm(ve) / norm(n));

        M = E - e * sin(E);
        a = 1 / (2 / norm(r) - norm(dr) * norm(dr) / m);

    };

    void get_all_Kepler() {
        std::cout << "KEPLER ELEMENTS\n" << "N " << N << "\n" << "i " << (i) << "\n" << "w " << w << "\n" << "a " << a
                  << "\n" << "e " << e << "\n" << "M " << M << "\n";
    };

    void get_all_Cartesian() {
        std::cout << "CARTESIAN COORDINATES\n" << "x " << x << "\n" << "y " << y << "\n" << "z " << z << "\n" << "vx "
                  << vx << "\n" << "vy " << vy << "\n" "vz " << vz << "\n";
    };

};



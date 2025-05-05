#include<iostream>
#include<cmath>
#include<vector>
#include<fstream>

double f(double x, double y){
    return (y/x) - 0.5 * sqrt(1 + pow(y/x, 2));
}

double runge_kutta(double x_0, double y_0, double h){
    double k_1 = f(x_0, y_0);
    double k_2 = f(x_0 + h/2, y_0 + h/2 * k_1);
    double k_3 = f(x_0 + h/2, y_0 + h/2 * k_2);
    double k_4 = f(x_0 + h, y_0 + h * k_3);

    return y_0 + (h / 6) * (k_1 + 2*k_2 + 2*k_3 + k_4);
}

int main(){
    double y_0 = 0;
    double h = -0.1;
    double y = 0;

    std::vector<std::vector<double>> trajectory;

    for (double x = 100; x >= 0; x += h) {
        y = runge_kutta(x, y_0, h);
        y_0 = y;
        trajectory.push_back({x, y});
    }

    std::ofstream outfile("trajectory.csv");
    if (!outfile.is_open()) {
        std::cerr << "Failed to open file for writing.\n";
        return 1;
    }

    outfile << "x,y\n";
    for (const auto& xy : trajectory) {
        outfile << xy[0] << "," << xy[1] << "\n";
    }

    outfile.close();
}

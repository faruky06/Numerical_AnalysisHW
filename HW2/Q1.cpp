#include<iostream>
#include<cmath>
#include<vector>
// Define function to check if a point is inside the circular disc
bool inside_circle(double x, double y) {
    return (pow(x - 0.25, 2) + pow(y - 0.25, 2)) < 0.125;
}

// Solve y numerically given x from the kidney equation
std::vector<double> solve_y(double x) {
    std::vector<double> y_values;
    for (double y = -0.3; y <= 1.0; y += 0.0001) { 
        double left = pow(x * x + y * y, 2);
        double right = pow(x, 3) + pow(y, 3);
        if (abs(left - right) < 1e-4) { // Numerical tolerance
            //if (!inside_circle(x, y)) { // Exclude points inside the disc
                y_values.push_back(y);
            //}
        }
    }
    return y_values;
}

// Rectangle method integration
double rectangle_method(double x_min, double x_max, double dx) {
    double area = 0.0;
    for (double x = x_min; x <= x_max; x += dx) {
        std::vector<double> y_vals = solve_y(x);
        if (y_vals.size() >= 2) { // Ensure we have upper and lower bounds
            double height = y_vals.back() + y_vals.front();
            area += height * dx;
        }
    }
    return area;
}

// Trapezoidal method integration
double trapezoidal_method(double x_min, double x_max, double dx) {
    double area = 0.0;  
    for (double x = x_min; x <= x_max - dx; x += dx) {
        std::vector<double> y_vals1 = solve_y(x);
        std::vector<double> y_vals2 = solve_y(x + dx);
        if (y_vals1.size() >= 2 && y_vals2.size() >= 2) {
            double height1 = y_vals1.back() + y_vals1.front();
            double height2 = y_vals2.back() + y_vals2.front();
            area += 0.5 * (height1 + height2) * dx;
        }
    }
    return area;
}

int main() {
    double x_min = -0.3, x_max = 1, dx = 0.0001; // Adjust range based on the graph

    double rect_area = rectangle_method(x_min, x_max, dx);
    double trap_area = trapezoidal_method(x_min, x_max, dx);

    std::cout << "Approximate area using Rectangle Method: " << rect_area << std::endl;
    std::cout << "Approximate area using Trapezoidal Method: " << trap_area << std::endl;

    return 0;
}
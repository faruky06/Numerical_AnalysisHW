#include <iostream>
#include <cmath>
#include <random>
#include <vector>
#include <fstream>
#define arrlen(a) (sizeof(a) / sizeof(a[0]))
int main(){
    double diam[13];
    //Arrays with specified diameter
    for(int i = 0; i < 10; i++){
        diam[i] = (i + 1) / 10.0;
    }
    diam[10] = 1.5; diam[11] = 2; diam[12] = 3;

    std::random_device rd;
    std::vector<std::vector<double>> probabilities;
    for(int i = 0; i < arrlen(diam); i++){
        double min_dist = diam[i] / 2; //minimum distance center of circle can be away from the edge
        std::uniform_real_distribution<double> dist(min_dist, 4 - min_dist); // generate a y-distance
        int touches[4] = {0,0,0,0}; //how many times touches (i) lines concurrently

        for(int k = 0; k < 4444444; k++){ 
            double y_value = dist(rd);
            int touch = 0; //how many lines circle touches in a single simulation
            for(int m = 1; m < 4; m++){ //checks distance of circle to line  y = 1, y = 2, y = 3
                if(abs(m - y_value) <= diam[i]/2)
                    touch++;
            }
            for(int m = 0; m < touch + 1; m++)
                touches[m]++;
        }
        std::vector<double> probabilities_i(4); 
        for(int k = 0; k < arrlen(touches); k++) //calculates probabilities
            probabilities_i[k] = static_cast<double>(touches[k]) / 4444444;
        probabilities.push_back(probabilities_i);
    }
    for (const auto& row : probabilities) {
        for (double val : row)
            std::cout << val << " ";
        std::cout << "\n";
    }
    
    std::ofstream outfile("circle_overlap.csv");
    if (!outfile.is_open()) {
    std::cerr << "Failed to open file for writing.\n";
    return 1;
    }

    // Optional: write header row
    outfile << "Diameter,P(0),P(1),P(2),P(3)\n";

    for (int i = 0; i < probabilities.size(); i++) {
        outfile << diam[i]; // First column: the diameter
        for (double val : probabilities[i]) {
        outfile << "," << val;
        }
        outfile << "\n";
    }

    outfile.close();
}
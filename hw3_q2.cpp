#include <iostream>
#include <cmath>
#include <algorithm>
#include <random>

const int SAMPLE_COUNT = 100000;
const double RECT_WIDTH = 1.0;
const double RECT_HEIGHT = 1.0 / sqrt(2);
const double PI = acos(-1);

//check if a point is inside the 4-leaf rose
bool isInsideRose(double px, double py) {
    double radius = sqrt(px * px + py * py);
    double angle = atan2(py, px);
    double roseRadius = fabs(sin(2 * angle));
    return radius <= roseRadius;
}
//check if a point is inside a rotated rectangle
bool isInsideRect(double px, double py, double cx, double cy, double rot) {
    double cosTheta = cos(rot);
    double sinTheta = sin(rot);
    double dx = px - cx;
    double dy = py - cy;
    double rotatedX = cosTheta * dx + sinTheta * dy;
    double rotatedY = -sinTheta * dx + cosTheta * dy;
    return (fabs(rotatedX) <= RECT_WIDTH / 2.0) && (fabs(rotatedY) <= RECT_HEIGHT / 2.0);
}
//Monte Carlo area estimation
double monteCarloArea(double cx, double cy, double rot, std::mt19937 &rng, std::uniform_real_distribution<> &squareDist) {
    int hits = 0;
    for (int i = 0; i < SAMPLE_COUNT; ++i) {
        double rx = squareDist(rng);
        double ry = squareDist(rng);
        if (isInsideRose(rx, ry) && isInsideRect(rx, ry, cx, cy, rot)) {
            hits++;
        }
    }
    return 4.0 * hits / SAMPLE_COUNT;
}
int main() {
    std::random_device entropySource;
    std::mt19937 shuffleEngine(entropySource());
    std::uniform_real_distribution<> coordDist(-1.0, 1.0);
    std::uniform_real_distribution<> angleDist(0.0, 2 * PI);
    double maxArea = 0.0;
    double optX = 0.0, optY = 0.0, optTheta = 0.0;
    for (int i = 0; i < 10000; ++i) {
        double testX = coordDist(shuffleEngine);
        double testY = coordDist(shuffleEngine);
        double testAngle = angleDist(shuffleEngine);
        double currentArea = monteCarloArea(testX, testY, testAngle, shuffleEngine, coordDist);
        if (currentArea > maxArea) {
            maxArea = currentArea;
            optX = testX;
            optY = testY;
            optTheta = testAngle;
        }
    }
    std::cout << "Best center: (" << optX << ", " << optY << ")\n";
    std::cout << "Best angle (radians): " << optTheta << "\n";
    std::cout << "Max estimated area: " << maxArea << "\n";
    return 0;
}
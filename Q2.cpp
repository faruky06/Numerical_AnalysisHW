#include<iostream>
#include<Eigen/Dense>

double func(const Eigen::VectorXd& coeff, int value){
    double eval;
    for(int i = 0; i < coeff.size(); i++){
        eval += coeff(i) * pow(value, coeff.size() - 1 - i);
    }
    return eval;
}
int main(){
    int rows = 5; int cols = 5;
    Eigen::MatrixXd x_values(rows, cols);
    Eigen::VectorXd y_values(5);
    for(int i = 0; i < rows; i++){
        for(int j = 0; j < cols; j++){
            x_values(i, j) = pow((rows - i), cols - 1 - j);
        }
    }
    y_values << 417, 398, 397, 407, 412; 
    Eigen::VectorXd solution_vec(5);
    solution_vec = x_values.colPivHouseholderQr().solve(y_values);;
    for(int i = 0; i < solution_vec.size(); i++){
        std::cout << solution_vec(i) << std::endl;
    }
}

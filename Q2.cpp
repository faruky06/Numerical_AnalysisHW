#include<iostream>
#include<Eigen/Dense>
//evaluates the function given the coefficient vector and the input to evaluate at
double func(const Eigen::VectorXd& coeff, int value){
    double eval;
    for(int i = 0; i < coeff.size(); i++){
        eval += coeff(i) * pow(value, coeff.size() - 1 - i);
    }
    return eval;
}
int main(){
    /*
    * Question 2.1
    * Initialize 5 x 5 matrix where each row has 
    */
    int rows = 5; int cols = 5;
    /* Matrix that hold target vectors x values evaluated to the (column_size - 1) power
    * Given by page 6 of topic 2 notes
    */
    Eigen::MatrixXd x_values(rows, cols); 
    Eigen::VectorXd target(5);
    for(int i = 0; i < rows; i++){
        for(int j = 0; j < cols; j++){
            x_values(i, j) = pow((rows - i), cols - 1 - j); 
        }
    }

    target << 417, 398, 397, 407, 412; //target vector with values to be interpolated
    Eigen::VectorXd solution_vec(5); //gives coefficients
    solution_vec = x_values.colPivHouseholderQr().solve(target); //Uses Eigen's QR decompisition and then solves to give solution vec
    std::cout<< "The interpolated polynomial is given by: ";
    for(int i = 0; i < solution_vec.size(); i++){
        if(i != solution_vec.size() - 1) 
            std::cout << solution_vec(i) << "x^" << (solution_vec.size() - 1 - i) << " + "; //outputs the polynomial to terminal
        else
            std::cout << solution_vec(i) << std::endl;
    }
    std::cout << "The interpolated polynomial evaluated at P(6) = " << func(solution_vec, 6) << std::endl;

    /*
    * Question 2.2
    * Quadratic Least Squares from given
    * 417       1 5 25      
    * 398       1 4 16      b
    * 397   =   1 3 9   *   m_1
    * 407       1 2 4       m_2
    * 412       1 1 1
    */
   Eigen::MatrixXd X(5, 3);
   for(int i = 0; i < X.rows(); i++){
    for(int j = 0; j < X.cols(); j++){
        X(i, j) = pow(cols - i, j); //populates matrix accordingly
    }
   }
    //Solves Beta^ = (X^T * X)(X^T) * y^ 
   Eigen::Vector3d quad_coeff = (X.transpose() * X).inverse() * X.transpose() * target; 
   std::cout << "The equation for a quadratic line of best fit: ";
   for(int i = quad_coeff.size() - 1; i >= 0; i--){
    if(i != 0)
        std::cout << quad_coeff(i) << "x^" << i << " + ";
    else
        std::cout << quad_coeff(i) << std::endl;
   }
   std::cout << "The quadratic fit evaluated at t = 6 is: " << func(quad_coeff.reverse(), 6);

}

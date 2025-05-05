#include<iostream>
#include<cmath>
#include<iomanip>
#include<random>
//#include<omp.h>
template<typename T1, typename T2, typename T3> //to encapsulate the root, iterations and fp operations
struct tuple{
    T1 t1;
    T2 t2;
    T3 t3;
    tuple(T1 t1, T2 t2, T3 t3) : t1(t1), t2(t2), t3(t3) {} 
    tuple(){}
};
struct df{ //struct to evaluate the derivative of the given function at a point x 
    long double operator() (long double x){
        return (-3 * pow(x, 2)) * (exp(-pow(x, 3))) - (4 * pow(x, 3)) - cos(x);
    }
};
struct f{
    df df;
    int iterations = 0; int floating_point_operations = 0;
   long double operator() (long double x){ //evaluates the function at a given point
        return exp(-pow(x, 3)) - pow(x, 4) - sin(x);
   }
   /*
   * takes interval to perform bisection on
   * returns tuple containing the root found within given tolerance, 
   * floating point operations and iterations done by method
   */
   tuple<long double, int, int> bisection(long double a, long double b){
    iterations = floating_point_operations = 0;
    tuple<long double, int, int> bisection_return; 
    long double c = ((a + b) / 2.0);
    while(!(abs(a-b) < 5e-6)){
        c = ((a + b) / 2.0);
        if((*this)(a) * (*this)(c) < 0)
            b = c;
        else
            a = c; 
        floating_point_operations += 3;
        iterations++;
    }
    bisection_return.t1 = ((a + b) / 2.0), bisection_return.t2 = floating_point_operations, bisection_return.t3 = iterations;
    return bisection_return;
   }
    /*
    *Newton's method given an initial guess, finds root within threshold
    */
    tuple<long double, int, int> newton(long double init_guess){
        iterations = 0; floating_point_operations = 3;
        tuple<long double, int, int> newton_return;
        long double x_n1 = init_guess - (((*this)(init_guess)) / df(init_guess)); //3 floating point operations
        while(!(abs((*this)(init_guess)) < (5e-6))){ //one evaluation of function is floating point operation
            init_guess = init_guess - (((*this)(init_guess)) / df(init_guess));
            floating_point_operations += 4;
            iterations++;
        }
        //assign values into tuple
        newton_return.t1 = init_guess; newton_return.t2 = floating_point_operations; newton_return.t3 = iterations; 
        return newton_return;
    }
    /*
    *Finds root within threshold given initial a_1, a_2
    */
    tuple<long double, int, int> secant(long double a, long double b){
        floating_point_operations = iterations = 0;
        tuple<long double, int, int> secant_return;
        long double new_b = 0;
        while(!(abs((*this)(a)) < 5e-6)){ //evaluation of function is a floating point operation
            new_b = a;
            a = a - (*this)(a) / (((*this)(a) - (*this)(b)) / (a - b)); // 8 operations evaluation of 3 functions, 2 divisions and 3 subtractions
            b = new_b; 
            floating_point_operations += 9; 
            iterations++;
        }
        secant_return.t1 = a; secant_return.t2 = floating_point_operations; secant_return.t3 = iterations;
        return secant_return;
    }
    /*
    *Randomly chooses values within given interval until one is within acceptable threshold of root 
    *Uses <random> instead of rand() function because they are more mathematically random 
    */
    tuple<long double, int, int> monte_carlo(long double lower, long double upper){
        floating_point_operations = iterations = 0;
        std::random_device rd;
        std::uniform_real_distribution<long double> dist(lower, upper);
        tuple<long double, int, int> mc_return;
        long double guess = dist(rd);
        while(!(abs((*this)(guess)) < 5e-6)){
            guess = dist(rd);
            iterations++;
            floating_point_operations += 2; // evaluation of function and generation of random number
        }
            mc_return.t1 = guess, mc_return.t2 = floating_point_operations, mc_return.t3 = iterations;
            return mc_return;
    }
    /*
    * Runs previous Monte Carlo algorithm with however many threads the program has access to
    * Stops when a thread finds a solution within the acceptable threshold
    */
   /*
    tuple<long double, int, int> monte_carlo_parallel(long double lower, long double upper){
        bool found_solution = false;
        floating_point_operations = iterations = 0;
        std::random_device rd;
        std::uniform_real_distribution<long double> dist(lower, upper);
        tuple<long double, int, int> result;
        #pragma omp parallel shared(found_solution) //found solution is shared between all threads
        {
        // local variables for each thread
        int iterations_thread = 0;  
        int floating_point_operations_thread = 0; 

        long double guess = dist(rd); //random guess
        while (!found_solution) { //if another thread has not found a solution execute code 
            guess = dist(rd);
            floating_point_operations_thread += 2;
            iterations_thread++;
            
            if (abs((*this)(guess)) < 5e-6) {
                #pragma omp atomic write //ensures only one thread can write to the memory of found_solution
                found_solution = true;
                result = tuple<long double, int, int>(guess, floating_point_operations_thread, iterations_thread);
            }
        }

        #pragma omp atomic 
        floating_point_operations += floating_point_operations_thread;
        
        #pragma omp atomic
        iterations += iterations_thread;
    }

    //return result with accumulated values
    result.t2 = floating_point_operations;  // total floating point operations
    result.t3 = iterations;  // total iterations
    return result;
    }
    */
};
int main(){
    f f;
    /*
    * Declares tuples containing all information for specificied initial conditions for each algo
    */
    tuple<long double, int, int> newton_return = f.newton(0);
    tuple<long double, int, int> secant_return = f.secant(-1, 1);
    tuple<long double, int, int> bisection_return = f.bisection(-1, 1);

    //keep track of runtime of parallel vs non parallel method, purely for personal use
    // has nothing to do with assignment, ignore parallel method.
    /*
    auto start_mcp = std::chrono::high_resolution_clock::now();
    tuple<long double, int, int> mc_return_parallel = f.monte_carlo_parallel(0.5, 0.75);
    auto end_mcp = std::chrono::high_resolution_clock::now();
    */

    auto start_mc = std::chrono::high_resolution_clock::now();
    tuple<long double, int, int> mc_return = f.monte_carlo(0.5, 0.75);
    auto end_mc = std::chrono::high_resolution_clock::now();

    std::chrono::duration<double> duration_mc = (end_mc - start_mc);
   //std::chrono::duration<double> duration_mcp = (end_mcp - start_mcp);
    long double time_per_op_mc = mc_return.t2 / duration_mc.count();
    //long double time_per_op_mcp = mc_return_parallel.t2 / duration_mcp.count();

    std::cout << "Bisection method with initial bracketing [-1, 1]: \n" << "        Root:"<< std::fixed << std::setprecision(8) << bisection_return.t1 << std::endl
    << "        floating point operations: " << bisection_return.t2 << std::endl
    << "        iterations: " << bisection_return.t3 << std::endl

    << "Secant method with initial points x0 = -1, x1 = 1: \n" << "     Root: " << secant_return.t1 << std::endl
    << "        floating point operations: " << secant_return.t2 << std::endl
    << "        iterations: " << secant_return.t3 << std::endl

    << "Newtons method with init guess x = 0: \n" << "      Root: " << newton_return.t1 << std::endl
    << "        floating point operations: " << newton_return.t2 << std::endl
    << "        iterations: " << newton_return.t3 << std::endl

    << "Monte Carlo method in range [0.5, 0.75]:\n" << "        Root: " << mc_return.t1 << std::endl
    << "        floating point operations: " << mc_return.t2 << std::endl
    << "        iterations: " << mc_return.t3 << std::endl
    << "        operations/second: " << time_per_op_mc << std::endl

    /*
    << "Monte Carlo parallel method in range [0.5, 0.75]:\n" << "       Root: " << mc_return_parallel.t1 << std::endl
    << "        floating point operations: " << mc_return_parallel.t2 << std::endl
    << "        iterations: " << mc_return_parallel.t3 << std::endl
    << "        operations/second: " << time_per_op_mcp;
    */
}

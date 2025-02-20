#include<iostream>
#include<cmath>
#include<iomanip>
#include<random>
#include<omp.h>
template<typename T1, typename T2, typename T3> //to encapsulate the root, iterations and fp operations
struct tuple{
    T1 t1;
    T2 t2;
    T3 t3;
    tuple(T1 t1, T2 t2, T3 t3) : t1(t1), t2(t2), t3(t3) {} // might not be best practice with naming  
    tuple(){}
};
struct df{
    long double operator() (long double x){
        return (-3 * pow(x, 2)) * (exp(-pow(x, 3))) - (4 * pow(x, 3)) - cos(x);
    }
};
struct f{
    df df;
    int iterations = 0; int floating_point_operations = 0;
   long double operator() (long double x){
        return exp(-pow(x, 3)) - pow(x, 4) - sin(x);
   }
   long double bisection(long double a, long double b){
    f f; // defines function from struct
    long double c = ((a + b) / 2.0);
    if(abs(a - b) < 1 * pow(10, -10))
        return c;
    else if (f(c)*f(a) < 0)
        return bisection(a, c);
    else
        return bisection(c, b);
    }

    tuple<long double, int, int> newton(long double init_guess){
        iterations = floating_point_operations = 0;
        tuple<long double, int, int> newton_return;
        long double x_n1 = init_guess - (((*this)(init_guess)) / df(init_guess));
        while(!(abs((*this)(init_guess)) < (1 * pow(10, -6)))){
            init_guess = init_guess - (((*this)(init_guess)) / df(init_guess));
            floating_point_operations++;
            iterations++;
        }
        //assign values into tuple
        newton_return.t1 = init_guess; newton_return.t2 = floating_point_operations; newton_return.t3 = iterations; 
        return newton_return;
    }
    long double secant(long double a, long double b){
        long double new_b = 0;
        long double new_a = 0;
        while(!((*this)(a) < 1e-6)){
            
        }
        return 0;
    }

    tuple<long double, int, int> monte_carlo(long double lower, long double upper){
        floating_point_operations = iterations = 0;
        std::random_device rd;
        std::uniform_real_distribution<long double> dist(lower, upper);
        tuple<long double, int, int> mc_return;
        long double guess = dist(rd);
        while(!(abs((*this)(guess)) < 1e-6)){
            guess = dist(rd);
            iterations++;
            floating_point_operations++;
        }
            mc_return.t1 = guess; mc_return.t2 = floating_point_operations; mc_return.t3 = iterations;
            return mc_return;
    }

    tuple<long double, int, int> monte_carlo_parallel(long double lower, long double upper){
        bool found_solution = false;
        floating_point_operations = iterations = 0;
        std::random_device rd;
        std::uniform_real_distribution<long double> dist(lower, upper);
        tuple<long double, int, int> result;
        #pragma omp parallel shared(found_solution)
        {
        int iterations_thread = 0;  // Local variable for each thread
        int floating_point_operations_thread = 0;  // Local variable for each thread

        long double guess = dist(rd);
        while (!found_solution) {
            guess = dist(rd);
            floating_point_operations_thread++;
            iterations_thread++;
            
            if (abs((*this)(guess)) < 1e-6) {
                #pragma omp atomic write
                found_solution = true;
                result = tuple<long double, int, int>(guess, floating_point_operations_thread, iterations_thread);
            }
        }

        #pragma omp atomic
        floating_point_operations += floating_point_operations_thread;
        
        #pragma omp atomic
        iterations += iterations_thread;
    }

    // Now return the result with accumulated values
    result.t2 = floating_point_operations;  // Total floating point operations
    result.t3 = iterations;  // Total iterations
    return result;
    }
};
int main(){
    f f;
    tuple<long double, int, int> newton_return = f.newton(0);
    auto start_mcp = std::chrono::high_resolution_clock::now();
    tuple<long double, int, int> mc_return_parallel = f.monte_carlo_parallel(-10, 5);
    auto end_mcp = std::chrono::high_resolution_clock::now();
    auto start_mc = std::chrono::high_resolution_clock::now();
    tuple<long double, int, int> mc_return = f.monte_carlo(-10, 5);
    auto end_mc = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration_mc = (end_mc - start_mc);
    std::chrono::duration<double> duration_mcp = (end_mcp - start_mcp);
    long double time_per_op_mc = mc_return.t2 / duration_mc.count();
    long double time_per_op_mcp = mc_return.t2 / duration_mcp.count();
    std::cout << "Bisection method with initial bracketing [-1, 1]: " << std::fixed << std::setprecision(16) << f.bisection(-1, 1) << std::endl
    << "Newtons method with init guess x = 0: \n" << "    Root: " << newton_return.t1 << std::endl
    << "    floating point operations: " << newton_return.t2 << std::endl

    << "    iterations: " << newton_return.t3 << std::endl
    << "Monte Carlo method in range [-10, 5]:\n" << "  Root: " << mc_return.t1 << std::endl
    << "    floating point operations: " << mc_return.t2 << std::endl
    << "    iterations: " << mc_return.t3 << std::endl
    << "    operations/second: " << time_per_op_mc << std::endl

    << "Monte Carlo parallel method in range [-10, 5]:\n" << "  Root: " << mc_return_parallel.t1 << std::endl
    << "    floating point operations: " << mc_return_parallel.t2 << std::endl
    << "    iterations: " << mc_return_parallel.t3 << std::endl
    << "    operations/second: " << time_per_op_mcp;
}
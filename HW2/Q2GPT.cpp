#include <iostream>
#include <random>
#include <cmath>

// Generic tuple structure
template <typename T1, typename T2, typename T3>
struct tuple {
    T1 t1;
    T2 t2;
    T3 t3;
    tuple(T1 t1, T2 t2, T3 t3) : t1(t1), t2(t2), t3(t3) {};
    tuple() {}
};

// Generic pair structure
template <typename T1, typename T2>
struct pair {
    T1 t1;
    T2 t2;
    pair(T1 t1, T2 t2) : t1(t1), t2(t2) {};
    pair() {}
};

// Function to allocate and populate an NxN matrix and a vector
tuple<double**, double*, int> populate_matrix_vec(int n) {
    double** matrix = new double*[n];
    double* vec = new double[n];

    std::uniform_real_distribution<double> dist(-1, 1);
    std::random_device rd;
    std::mt19937 gen(rd());

    for (int i = 0; i < n; i++) {
        vec[i] = 1;  // Initialize solution vector to [1,1,1,...,1]
        matrix[i] = new double[n];
        for (int j = 0; j < n; j++) {
            matrix[i][j] = dist(gen);  // Populate matrix with random values
        }
    }

    return tuple<double**, double*, int>(matrix, vec, n);
}

// Function to create a copy of a matrix
double** copy_matrix(double** original, int n) {
    double** copy = new double*[n];
    for (int i = 0; i < n; i++) {
        copy[i] = new double[n];
        for (int j = 0; j < n; j++) {
            copy[i][j] = original[i][j];
        }
    }
    return copy;
}

// Function to safely delete a dynamically allocated matrix
void delete_matrix(double** matrix, int n) {
    if (!matrix) return;
    for (int i = 0; i < n; i++) {
        delete[] matrix[i];
    }
    delete[] matrix;
}

// Function to perform Gaussian elimination and return RREF matrix + solution vector
pair<double**, double*> gaussian_elim_X_unknown(tuple<double**, double*, int>& input) {
    int n = input.t3;

    // Make a copy of the matrix so the original is not modified
    double** mat = copy_matrix(input.t1, n);
    double* vec = new double[n];

    // Copy vector to prevent modification
    for (int i = 0; i < n; i++) {
        vec[i] = input.t2[i];
    }

    for (int i = 0; i < n; i++) {
        // Find the pivot row
        int maxRow = i;
        for (int k = i + 1; k < n; k++) {
            if (fabs(mat[k][i]) > fabs(mat[maxRow][i])) {
                maxRow = k;
            }
        }
        // Swap rows if needed
        if (i != maxRow) {
            std::swap(mat[i], mat[maxRow]);
            std::swap(vec[i], vec[maxRow]);
        }

        // Scale row to make pivot = 1
        double pivot = mat[i][i];
        if (fabs(pivot) < 1e-10) continue;
        for (int j = 0; j < n; j++) {
            mat[i][j] /= pivot;
        }
        vec[i] /= pivot;

        // Make other entries in the column zero
        for (int k = 0; k < n; k++) {
            if (k == i) continue;
            double factor = mat[k][i];
            for (int j = 0; j < n; j++) {
                mat[k][j] -= factor * mat[i][j];
            }
            vec[k] -= factor * vec[i];
        }
    }
    return pair<double**, double*>(mat, vec);
}

// Function to compute matrix-vector multiplication
double* dot_product(double** matrix, double* vec, int n) {
    double* result = new double[n];

    for (int i = 0; i < n; i++) {
        result[i] = 0;  // Initialize sum to 0
        for (int j = 0; j < n; j++) {
            result[i] += matrix[i][j] * vec[j];
        }
    }

    return result;
}

int main() {
    int n = 66;  // Matrix size

    // Generate a random matrix and solution vector
    tuple<double**, double*, int> mat_vec = populate_matrix_vec(66);

    // Print original matrix
    std::cout << "Original Matrix:\n";
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            std::cout << mat_vec.t1[i][j] << " ";
        }
        std::cout << std::endl;
    }

    // Perform Gaussian elimination (pass by reference)
    pair<double**, double*> reduced = gaussian_elim_X_unknown(mat_vec);
    double** mat_rref = reduced.t1;
    double* sol_vec = reduced.t2;

    // Print RREF matrix
    std::cout << "\nRREF Matrix:\n";
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            std::cout << mat_rref[i][j] << " ";
        }
        std::cout << std::endl;
    }

    // Print solution vector
    std::cout << "\nSolution Vector (should be [1,1,1,1,1]):\n";
    for (int i = 0; i < n; i++) {
        std::cout << sol_vec[i] << " ";
    }
    std::cout << std::endl;

    // Verify A * x = b
    double* check = dot_product(mat_vec.t1, sol_vec, n);

    std::cout << "\nChecking A * x:\n";
    for (int i = 0; i < n; i++) {
        std::cout << check[i] << " ";
    }
    std::cout << std::endl;

    // Clean up allocated memory
    delete_matrix(mat_vec.t1, n);
    delete_matrix(mat_rref, n);
    delete[] mat_vec.t2;
    delete[] sol_vec;
    delete[] check;

    return 0;
}
#include <iostream>
#include <vector>
#include <math.h>
#include <chrono>

using namespace std;

float det(vector<vector<float>> A) {
    if (A.size() == 2) {
        return (A[0][0] * A[1][1]) - (A[0][1] * A[1][0]);
    }
    float d = 0;
    int c;
    for (int i=0 ; i < A.size() ; i++) {
        vector<vector<float>> B((A.size())-1 , vector<float>((A.size())-1));
        for (int j=1 ; j < A.size() ; j++) {
            c = 0;
            for (int k=0 ; k < A.size() ; k++) {
                if (k != i) {
                    B[j-1][c] = A[j][k];
                    c++;
                }
            }
        }
        d += pow(-1 , i) * A[0][i] * det(B);
    }
    return d;
}

vector<vector<float>> recursive_inverse(vector<vector<float> > A) {
    int m,n;
    float d = det(A);

    vector<vector<float>> B(A.size() , vector<float>(A[0].size()));
    for (int i = 0 ; i < A.size() ; i++) {
        for (int j = 0 ; j < A.size() ; j++) {
            vector<vector<float>> C((A.size())-1 , vector<float>((A[0].size())-1));
            m = 0;
            n = 0;

            for (int r=0 ; r < A.size() ; r++) {
                if (r != i) {
                    n=0;
                    for (int c = 0 ; c < A.size() ; c++) {
                        if (c != j) {
                            C[m][n] = A[r][c];
                            n++;
                        }

                    }
                    m++;
                }
            }
            B[i][j] = pow (-1 , i+j) * det(C);
        }
    }

    vector<vector<float>> D(A.size() , vector<float>(A[0].size()));
    for (int i=0 ; i < A.size() ; i++) {
        for (int j=0 ; j < A.size() ; j++) {
            D[i][j] = B[j][i];
        }
    }

    vector<vector<float>> E(A.size() , vector<float>(A[0].size()));
    for (int i=0 ; i < A.size() ; i++) {
        for (int j=0 ; j < A.size() ; j++) {
            E[i][j] = D[i][j] / d;
        }
    }

    return E;
}

vector<vector<float> > gauss_inv(vector<vector<float> > A, double epsilon = 1e-5) {
    vector<vector<float> > B(A.size(), vector<float>(2 * A[0].size(), 0));
    vector<vector<float> > C(A.size(), vector<float>(A[0].size(), 0));
    float d, f;
    for (int i = 0; i < A.size(); i++) {
        for (int j = 0; j < A[0].size(); j++) {
            B[i][j] = A[i][j];
        }
        B[i][i + (A[0].size())] = 1;
    }
    for (int i = 0; i < A.size(); i++) {
        d = B[i][i];
        for (int j = 0; j < 2 * A[0].size(); j++) {
            B[i][j] = B[i][j] / d;
        }
        for (int k = 0; k < A.size(); k++) {
            if (k != i) {
                f = B[k][i];
                for (int j = 0; j < 2 * (A[0].size()); j++) {
                    B[k][j] = B[k][j] - (f * B[i][j]);
                }
            }
        }
    }
    for (int i = 0; i < A.size(); i++) {
        for (int j = 0; j < A.size(); j++) {
            C[i][j] = B[i][j + (A.size())];
        }
    }
    return C;
}

vector<vector<float> > multi(vector<vector<float> > A, vector<vector<float> > B, double epsilon = 1e-5) {
    vector<vector<float> > R(A.size(), vector<float>(B[0].size(), 0));
    for (int i = 0; i < A.size(); i++) {
        for (int j = 0; j < B[0].size(); j++) {
            for (int k = 0; k < B.size(); k++) {
                R[i][j] += A[i][k] * B[k][j];
            }
            if (R[i][j] < epsilon) {
                R[i][j] = 0;
            } else if (1 - epsilon < R[i][j] < 1 + epsilon) {
                R[i][j] = 1;
            }
        }
    }
    return R;
}

int main() {
    vector<vector<float>> A = {{3, 1, 4, 1, 5, 9, 2, 6, 5, 3},
                               {5, 8, 9, 7, 9, 3, 2, 3, 8, 4},
                               {6, 2, 6, 4, 3, 3, 8, 3, 2, 7},
                               {9, 5, 0, 2, 8, 8, 4, 1, 9, 7},
                               {1, 6, 9, 3, 9, 7, 5, 1, 0, 5},
                               {8, 2, 6, 5, 3, 5, 8, 9, 7, 9},
                               {3, 2, 3, 8, 4, 6, 2, 6, 4, 3},
                               {3, 8, 3, 2, 7, 9, 5, 0, 2, 8},
                               {8, 4, 1, 9, 7, 1, 6, 9, 3, 9},
                               {7, 5, 1, 0, 5, 8, 2, 6, 5, 3}};

    auto recursive_start = chrono::high_resolution_clock::now();
    vector<vector<float>> RB = recursive_inverse(A);
    vector<vector<float>> RC = multi(A,RB);

    cout << "Recursive Inverse:" << endl;
    for (int i = 0 ; i < A.size() ; i++) {
        for (int j = 0 ; j < A.size() ; j++) {
            cout << RB[i][j] << " ";
        }
        cout << endl;
    }



    cout << "*************************************************************************" << endl;

    auto  recursive_end = chrono::high_resolution_clock::now();
    auto diff = chrono::duration_cast<chrono::microseconds>(recursive_end - recursive_start);
    cout << " Recursive inverse runtime: " << diff.count() << "microseconds" << endl;

    cout << "*************************************************************************************" << endl;

    auto gauss_start = chrono::high_resolution_clock::now();
    vector<vector<float> > GB = gauss_inv(A);
    vector<vector<float> > GC = multi(A, GB);



    cout << "gauss jordan:" << endl;
    for (int i = 0; i < A.size(); i++) {
        for (int j = 0; j < A.size(); j++) {
            cout << GB[i][j] << " ";
        }
        cout << endl;
    }

    cout << "*************************************************************************" << endl;

    auto gauss_end = chrono::high_resolution_clock::now();
    auto result = chrono::duration_cast<chrono::microseconds>(gauss_end-gauss_start);

    cout << "Gauss jordan runtime: " << result.count() << "microseconds" << endl;

    return 0;
}
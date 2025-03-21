#include <iostream>
#include <vector>
#include <limits>
#include <algorithm>
#include <iomanip>

using namespace std;

void printSimplexTable(const vector<vector<double>>& A, const vector<double>& b, const vector<double>& c, double v, const vector<int>& N, const vector<int>& B) {
    int m = A.size();
    int n = A[0].size();

    cout << "Базис\t";
    for (int j = 0; j < n; ++j) {
        cout << "x" << N[j] << "\t";
    }
    cout << "Свободный\n";

    for (int i = 0; i < m; ++i) {
        cout << "x" << B[i] << "\t";
        for (int j = 0; j < n; ++j) {
            cout << fixed << setprecision(6) << A[i][j] << "\t";
        }
        cout << fixed << setprecision(6) << b[i] << "\n";
    }

    cout << "Целевая функция\t";
    for (int j = 0; j < n; ++j) {
        cout << fixed << setprecision(6) << c[j] << "\t";
    }
    cout << fixed << setprecision(6) << v << "\n\n";
}
//l - строка для поворота
//vhodperemen - столбец для поворота
void change(vector<int>& N, vector<int>& B, vector<vector<double>>& A, vector<double>& b, vector<double>& c, double& v, int l, int vhodperemen) {   //ПОВОРОТ
    int m = A.size();
    int n = A[0].size();

    b[l] /= A[l][vhodperemen];  //A[l][vhodperemen] - ведущий элемент
    for (int j = 0; j < n; ++j) {  
        if (j != vhodperemen) A[l][j] /= A[l][vhodperemen];  //A[l][vhodperemen] должен стать 1
    }
    A[l][vhodperemen] = 1 / A[l][vhodperemen];

    for (int i = 0; i < m; ++i) {  //обновление значений для всех строк кроме L
        if (i != l) {
            b[i] -= A[i][vhodperemen] * b[l];
            for (int j = 0; j < n; ++j) {
                if (j != vhodperemen) A[i][j] -= A[i][vhodperemen] * A[l][j];
            }
            A[i][vhodperemen] = -A[i][vhodperemen] * A[l][vhodperemen];
        }
    }

    v += c[vhodperemen] * b[l];     //обновление коэффициентов целевой функции
    for (int j = 0; j < n; ++j) {
        if (j != vhodperemen) c[j] -= c[vhodperemen] * A[l][j];
    }
    c[vhodperemen] = -c[vhodperemen] * A[l][vhodperemen];

    swap(N[vhodperemen], B[l]);  //Меняем местами небазисную переменную N[vhodperemen] и базисную переменную B[l]. 
                     //Это отражает тот факт, что переменная vhodperemen теперь входит в базис, а переменная l выходит из базиса
}


void simplex(vector<vector<double>> A, vector<double> b, vector<double> c, bool isMaximization) {
    try {
        int m = A.size();
        int n = A[0].size();
        vector<int> N(n); 
        vector<int> B(m);
        double v = 0; 

        for (int i = 0; i < n; ++i) N[i] = i + 1;
        for (int i = 0; i < m; ++i) B[i] = n + i + 1;

        cout << "Исходная таблица:\n";
        printSimplexTable(A, b, c, v, N, B);

        while (true) {
            int vhodperemen = -1;   //для улучшения целевой функции
            for (int j = 0; j < n; ++j) {
                if (c[j] > 0) {
                    vhodperemen = j;
                    break;
                }
            }
            if (vhodperemen == -1) break;

            vector<double> delta(m, numeric_limits<double>::infinity());   //выбор исключаемой переменной
            for (int i = 0; i < m; ++i) {
                if (A[i][vhodperemen] > 0) {
                    delta[i] = b[i] / A[i][vhodperemen];
                }
            }

            int l = distance(delta.begin(), min_element(delta.begin(), delta.end()));
            if (delta[l] == numeric_limits<double>::infinity()) {
                return;
            }
            change(N, B, A, b, c, v, l, vhodperemen);

            
        }

        cout << "Оптимальное решение:\n";
        for (int i = 0; i < n + m; ++i) {
            if (find(B.begin(), B.end(), i + 1) != B.end()) {
                int index = distance(B.begin(), find(B.begin(), B.end(), i + 1));
                cout << "x" << i + 1 << ": " << b[index] << endl;
            }
            else {
                cout << "x" << i + 1 << ": " << 0 << endl;
            }
        }
        cout << "Оптимальное значение: " << (isMaximization ? v : -v) << endl;
    }
    catch (const char* msg) {
        cout << msg << endl;
    }
}


int main() {
    setlocale(LC_ALL, "RUS");
    vector<vector<double>> A = {
        {1.01, 1.01, 9.45, 27},
        {0.002, 0.00166, 0, 0.0018},
        {0, 0, 0.033, 0},
        {0, 1, 0, 1}
    };
    vector<double> b = { 14000, 21, 16, 800 };
    vector<double> c = { 2.4, 2.7, 13.8, 0.77 };
    simplex(A, b, c, true);

    return 0;
}
#include <fstream>
#include <iostream>
#include <vector>
#include <cmath>
using namespace std;

extern "C"
{
    void dsolve_(int* neq, double* rpar, int* ipar, double* y, double* t, double* tstep);
}

//проверка на корректность ввода (на каждой стадии одинаковое количество элементов)
//fix - фиксированное количество (количество элементов на первой стадии)
//current - количество на текущей стадии
bool size_check(int fix, int current);

//переводим двумерный массив в одномерный для передачи в функцию f
double* to_rpar(vector<vector<double>> vector_2d);

int main()
{
    bool is_correct = true;
    //храним фиксированный размер (до первой стадии равен -1)
    int fix_size = -1;

    //файловый поток для ввода начальных данных
    ifstream fin;
    fin.open("reactions.txt", ios::in);

    //массив, в котором будем хранить стехиометрические коэффиценты
    vector<vector<double>> reaction;
    int i = 1;
    //флаг следит за тем, что мы доходим до конца строки
    bool flag = false;
    vector<double> temp;
    while (!fin.eof())
    {
        while ((fin.peek() != '\n' && !fin.eof()) || flag)
        {
            double coeff;
            fin >> coeff;
            temp.push_back(coeff);
            if (fin.peek() != '\n') flag = false;
        }

        //сейчас указатель в файле стоит на конце строки, поэтому для дальнейшего чтения в цикле выше нужно отметить, что сейчас начинается запись новой строки
        flag = true;

        //проверка, что не было ошибок в вводе
        if (fix_size < 0) fix_size = temp.size();
        else if (!size_check(fix_size, temp.size()))
        {
            cout << "Ошибка ввода, строка" << i << ": проверьте корректность начальных данных" << endl;
            is_correct = false;
            break;
        }

        reaction.push_back(temp);
        temp.clear();
        bool eof = fin.eof();
        i++;
    }

    if (is_correct)
    {
        int n = (fix_size - 1) / 2;
        double* c = new double[n];
        for (int i = 0; i < n; i++) c[i] = i == 0 ? 1. : 0.;
        double* rpar = to_rpar(reaction);
        int* ipar = new int[1];
        ipar[0] = reaction.size();
        double t0 = 0., dt = 0.01;
        dsolve_(&n, rpar, ipar, c, &t0, &dt);
    }
}

bool size_check(int fix, int current)
{
    return fix == current;
}

double* to_rpar(vector<vector<double>> vector_2d)
{
    int N = vector_2d.size();
    int M2 = vector_2d[0].size();
    double* rpar = new double[N * M2];
    for (int n = 0; n < N; n++)
    {
        for (int i = 0; i < M2; i++)
            {
                rpar[n * M2 + i] = vector_2d[n][i];
            }
    }

    return rpar;
}

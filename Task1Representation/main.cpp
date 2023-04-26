// Erokhin Evgenii DSAI-03
#include <iostream>
#include <vector>
#include <iomanip>
#include <cmath>
#include "cstdio";
#define GNUPLOT_NAME "C:\\gnuplot\\bin\\gnuplot -persist"



using namespace std;
class Matrix{
public:


    vector<vector<double>> matrixNew;
    int row;
    int column;
    Matrix(){
        row = 0;
        column = 0;
    }
    Matrix(int r, int c) {
        row = r;
        column = c;
        matrixNew.resize(row,vector<double>(column));

    }

    void add(int i, int j, double value)  { //add element
        matrixNew[i][j] = value;
    }

    double get(int i, int j)  { //get element
        return matrixNew[i][j];
    }

    Matrix operator=( Matrix& M) {
        matrixNew=M.matrixNew;
        row = M.row;
        column = M.column;
        return *this;
    }

    Matrix operator+(Matrix& M) {
        Matrix newM;
        try {
            if (row != M.row || column != M.column) {
                exception e;
                throw exception(e);
            } else{
                newM = *this;
                for (int i = 0; i < row; i++)
                    for (int j = 0; j < column; j++)
                        newM.add(i, j, matrixNew[i][j] + M.get(i, j));  //overload +

            }
        } catch (exception e){cout << "Error: the dimensional problem occurred\n";}
        return newM;
    }

    Matrix operator-(Matrix& M)  {Matrix newM;
        try {
            if (row != M.row || column != M.column) {
                exception e;
                throw exception(e);
            } else{ newM = *this;
                for (int i = 0; i < row; i++)
                    for (int j = 0; j < column; j++)
                        newM.add(i, j, matrixNew[i][j] - M.get(i, j)); //overload -

            }
        } catch (exception e){cout << "Error: the dimensional problem occurred\n";}
        return newM;
    }

    // overload of operator "*"
    Matrix operator*(Matrix& M)  {
        Matrix forSwitch;
        try {
            if (column != M.row) {
                exception e;
                throw exception(e);
            } else{ Matrix temp = Matrix(row, M.column);
                forSwitch = temp;
                for (int i = 0; i < row; i++) {
                    for (int j = 0; j < M.column; j++) {  //overload *
                        forSwitch.add(i, j, 0);
                        for (int k = 0; k < M.row; k++)
                            forSwitch.add(i, j, forSwitch.get(i, j)
                                                + matrixNew[i][k] * M.get(k, j));
                    }
                }
            }
        } catch (exception e){cout << "Error: the dimensional problem occurred\n";}
        return forSwitch;
    }

    Matrix transpose()  {
        Matrix forSwitch(column, row);
        for(int i = 0; i < row; i++)
            for(int j = 0; j < column; j++)
                forSwitch.add(j, i, matrixNew[i][j]); //transpose matrix
        return forSwitch;
    }
};
class SquareMatrix: public Matrix {
public:



    SquareMatrix() {
        row = 0;
        column = 0;
    }

    explicit SquareMatrix(int r) {
        row = column = r;
        matrixNew.resize(row, vector<double>(column));

    }

    SquareMatrix operator=(SquareMatrix &M) {
        Matrix *m1 = &M;
        Matrix *sq = this;
        Matrix res = sq->operator=(*m1);
        SquareMatrix *res1 = (SquareMatrix *) &res;
        return *res1;

    }

    SquareMatrix operator*(SquareMatrix &M) {
        Matrix *m1 = &M;
        Matrix *sq = this;
        Matrix res = sq->operator*(*m1);
        SquareMatrix *res1 = (SquareMatrix *) &res;
        return *res1;
    }

    SquareMatrix operator-(SquareMatrix &M) {
        Matrix *m1 = &M;
        Matrix *sq = this;
        Matrix res = sq->operator-(*m1);
        SquareMatrix *res1 = (SquareMatrix *) &res;
        return *res1;
    }

    SquareMatrix operator+(SquareMatrix &M) {
        Matrix *m1 = &M;
        Matrix *sq = this;
        Matrix res = sq->operator+(*m1);
        SquareMatrix *res1 = (SquareMatrix *) &res;
        return *res1;
    }
};
class IdentityMatrix: public SquareMatrix {

public:

    IdentityMatrix() {
        row = 0;
        column = 0;
    }

    explicit IdentityMatrix(int r) {
        row = column = r;
        matrixNew.resize(row, vector<double>(column));
        for (int ro = 0; ro < row; ro++) {
            for (int col = 0; col < column; col++) {
                if (ro == col) {
                    matrixNew[ro][col] = 1;
                } else {
                    matrixNew[ro][col] = 0;
                }
            }
        }

    }

    IdentityMatrix operator=(IdentityMatrix &M) {
        SquareMatrix *m1 = &M;
        SquareMatrix *sq = this;
        SquareMatrix res = sq->operator=(*m1);
        IdentityMatrix *res1 = (IdentityMatrix *) &res;
        return *res1;

    }

    IdentityMatrix operator*(IdentityMatrix &M) {
        SquareMatrix *m1 = &M;
        SquareMatrix *sq = this;
        SquareMatrix res = sq->operator*(*m1);
        IdentityMatrix *res1 = (IdentityMatrix *) &res;
        return *res1;
    }

    IdentityMatrix operator-(IdentityMatrix &M) {
        SquareMatrix *m1 = &M;
        SquareMatrix *sq = this;
        SquareMatrix res = sq->operator-(*m1);
        IdentityMatrix *res1 = (IdentityMatrix *) &res;
        return *res1;
    }

    IdentityMatrix operator+(IdentityMatrix &M) {
        SquareMatrix *m1 = &M;
        SquareMatrix *sq = this;
        SquareMatrix res = sq->operator+(*m1);
        IdentityMatrix *res1 = (IdentityMatrix *) &res;
        return *res1;
    }

};
class EliminationMatrix:public IdentityMatrix {
public:
    EliminationMatrix(){
        row = 0;
        column = 0;
    }


    explicit EliminationMatrix(SquareMatrix& A, int i,int j) {
        IdentityMatrix temp(A.row);
        matrixNew.resize(A.row, vector<double>(A.column));

        SquareMatrix Anew = A;
        int ro = 0;
        double element;
        element = Anew.matrixNew[i][j]/Anew.matrixNew[j][j];
        temp.matrixNew[i][j] = -1*element;
        EliminationMatrix *res1 = (EliminationMatrix *) &temp;
        *this = *res1;

    }


};
class PermutationMatrix:public IdentityMatrix {
public:
    PermutationMatrix(){
        column = 0;
        row = 0;
    }


    explicit PermutationMatrix(SquareMatrix& A, int i,int j) {
        IdentityMatrix temp(A.row);
        matrixNew.resize(A.row, vector<double>(A.column));

        SquareMatrix Anew = A;
        vector<double> temp1 = temp.matrixNew[i];
        temp.matrixNew[i] = temp.matrixNew[j];
        temp.matrixNew[j] = temp1;
        PermutationMatrix *res1 = (PermutationMatrix *) &temp;
        *this = *res1;

    }


};

static SquareMatrix Inverse(SquareMatrix&A){
    IdentityMatrix temp(A.row);
    Matrix ansToPrint(A.row,A.column+A.column);
    for (int i = 0; i < A.row; i++) {
        for (int j = 0; j < A.column; j++) {
            ansToPrint.matrixNew[i][j] = A.matrixNew[i][j];
        }
    }
    for (int i = 0; i < A.row; i++) {
        for (int j = A.column; j < ansToPrint.column; j++) {
            ansToPrint.matrixNew[i][j] = temp.matrixNew[i][j-A.column];
        }
    }
    for (int j = 0; j < A.column; j++) {

        int fabsMax = -1111111;
        int rowMax = -1;
        for (int i = j; i < A.row; i++) {
            if (fabs(A.matrixNew[i][j]) > fabsMax) {
                fabsMax = fabs(A.matrixNew[i][j]);
                rowMax = i;
            }
        }
        if (j != rowMax) {
            PermutationMatrix x(A,j,rowMax);
            IdentityMatrix V = x;
            SquareMatrix Z = x;
            SquareMatrix res1 = Z * A;
            IdentityMatrix res2 = V * temp;
            temp = res2;
            A = res1;
            for (int i = 0; i < A.row; i++) {
                for (int j = 0; j < A.column; j++) {
                    ansToPrint.matrixNew[i][j] = A.matrixNew[i][j];
                }
            }
            for (int i = 0; i < A.row; i++) {
                for (int j = A.column; j < ansToPrint.column; j++) {
                    ansToPrint.matrixNew[i][j] = temp.matrixNew[i][j-A.column];
                }
            }
        }
        for (int k = j+1;k<A.row;k++){
            if (A.matrixNew[k][j]!=0) {
                EliminationMatrix H(A, k, j);
                IdentityMatrix V = H;
                SquareMatrix Y = H;
                SquareMatrix res = Y * A;
                IdentityMatrix res1 = V*temp;
                A = res;
                temp = res1;
                for (int i = 0; i < A.row; i++) {
                    for (int j = 0; j < A.column; j++) {
                        ansToPrint.matrixNew[i][j] = A.matrixNew[i][j];
                    }
                }
                for (int i = 0; i < A.row; i++) {
                    for (int j = A.column; j < ansToPrint.column; j++) {
                        ansToPrint.matrixNew[i][j] = temp.matrixNew[i][j-A.column];
                    }
                }
            }
        }
    }

    for (int j = A.column-1; j >= 0; j--) {
        for (int k = j - 1; k >=0; k--) {
            if (A.matrixNew[k][j] != 0) {
                EliminationMatrix H(A, k, j);
                IdentityMatrix V = H;
                SquareMatrix Y = H;
                SquareMatrix res = Y * A;
                IdentityMatrix res1 = V*temp;
                A = res;
                temp = res1;
                for (int i = 0; i < A.row; i++) {
                    for (int j = 0; j < A.column; j++) {
                        ansToPrint.matrixNew[i][j] = A.matrixNew[i][j];
                    }
                }
                for (int i = 0; i < A.row; i++) {
                    for (int j = A.column; j < ansToPrint.column; j++) {
                        ansToPrint.matrixNew[i][j] = temp.matrixNew[i][j-A.column];
                    }
                }
            }
        }
    }
    int pos1=0,pos2=0;
    for (int i = 0; i < A.row; i++) {
        for (int j = 0; j < A.column; j++) {
            temp.matrixNew[i][j]=temp.matrixNew[i][j]/A.matrixNew[pos1][pos2];
        }
        A.matrixNew[pos1][pos2]=1;
        pos1++;
        pos2++;
    }

    SquareMatrix answer = temp;
    for (int i = 0; i < A.row; i++) {
        for (int j = 0; j < A.column; j++) {
            ansToPrint.matrixNew[i][j] = A.matrixNew[i][j];
        }
    }
    for (int i = 0; i < A.row; i++) {
        for (int j = A.column; j < ansToPrint.column; j++) {
            ansToPrint.matrixNew[i][j] = temp.matrixNew[i][j-A.column];
        }
    }
    return answer;
};

istream& operator>>(istream& in, Matrix& matrix) {
    double val;
    for (int i = 0; i < matrix.row; i++) {
        for (int j = 0; j < matrix.column; j++) {  //overload >>
            in >> val;
            matrix.add(i, j, val);
        }
    }
    return in;
}
ostream& operator<<(ostream& out, Matrix& matrix) {
    for (int i = 0; i < matrix.row; i++) {
        for (int j = 0; j < matrix.column; j++) {  //overload <<
            out <<fixed << setprecision(4)<<matrix.get(i, j) << " ";
        }
        out << endl;
    }
    return out;
}


int main(int argc,  char * argv[]) {

    FILE* pipe = _popen(GNUPLOT_NAME, "w");



    int lengthOfDataSet;
    cin >> lengthOfDataSet;   // read length of data set
    vector<double> tI;  //vector of t(i)
    vector<double> bI; //vector of b(i)
    for (int i=0;i<lengthOfDataSet;i++){  //fill vectors
        double temp1;
        cin>>temp1;
        tI.push_back(temp1);
        double temp2;
        cin>>temp2;
        bI.push_back(temp2);
    }
    int degree;
    cin >> degree;
    degree = degree+1;
    Matrix A(lengthOfDataSet,degree);
    for (int i=0;i<lengthOfDataSet;i++){
        for (int j=0;j<degree;j++){
            A.matrixNew[i][j] = pow(tI[i],j); //fill matrix A
        }
    }

    Matrix A_TA = A.transpose()*A; // A(transpose)*A


    // (A(transpose)*A) Inverse
    SquareMatrix *res1 = (SquareMatrix *) &A_TA;
    SquareMatrix A_TAInverse = Inverse(*res1);

    // A(transpose)*A * b
    Matrix b(lengthOfDataSet,1);
    for (int i=0;i<lengthOfDataSet;i++){
        b.matrixNew[i][0] = bI[i];
    }
    Matrix A_Tb = A.transpose()*b;


    Matrix res2 = A_TAInverse;
    Matrix x = res2*A_Tb;

    fprintf(pipe, "plot [-30 : 30] [-30 : 30] %lf*x**3 + %lf*x**2 + %lf*x**1 + %lf*x**0 , '-' using 1:2 with points\n",
            x.matrixNew[3][0], x.matrixNew[2][0], x.matrixNew[1][0], x.matrixNew[0][0]);
    fprintf(pipe, "plot '-' w p ls 1\n");
    for (int i = 0; i < lengthOfDataSet; i++) {
        fprintf(pipe, "%f\t%f\n", tI[i], bI[i]);
    }


}



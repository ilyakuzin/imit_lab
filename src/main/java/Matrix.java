import java.io.FileOutputStream;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Arrays;

public class Matrix {

    private double[][] matrix;
    private int maxLineElemNumber;
    private int maxColumnElemNumber;
    private double maxElem;
    private double[][] inverseMatrix;
    private double norm;

    public Matrix(double[][] matrix) {
        this.matrix = matrix;
    }

    //поиск обратной матрицы
    public double[][] getInverseMatrix() {

        double[][] mCopy = new double[matrix.length][matrix.length];
        double[][] matrixToChange = new double[matrix.length][matrix.length];
        double me;

        //сохраняем матрицу
        for (int i = 0; i < getMatrix().length; i++) {
            for (int j = 0; j < getMatrix().length; j++) {
                mCopy[i][j] = matrix[i][j];
            }
        }

        //задаю единичную матрицу
        for (int i = 0; i < getMatrix().length; i++) {
            for (int j = 0; j < getMatrix().length; j++) {
                if (i == j) {
                    matrixToChange[i][j] = 1;
                } else {
                    matrixToChange[i][j] = 0;
                }
            }
        }

        //делим на главный коэффицент
        for (int k = 0; k < getMatrix().length; k++) {
            me = matrix[k][k];
            for (int j = 0; j < getMatrix().length; j++) {
                matrix[k][j] /= me;
                matrixToChange[k][j] /= me;
            }

            //прямой ход
            for (int i = k + 1; i < getMatrix().length; i++) {
                me = matrix[i][k];
                for (int j = 0; j < getMatrix().length; j++) {
                    matrix[i][j] -= matrix[k][j] * me;
                    matrixToChange[i][j] -= matrixToChange[k][j] * me;
                }
            }
        }
        //обратный
        for (int k = getMatrix().length - 1; k > 0; k--) {
            for (int i = k - 1; i >= 0; i--) {
                me = matrix[i][k];
                for (int j = 0; j < getMatrix().length; j++) {
                    matrix[i][j] -= matrix[k][j] * me;
                    matrixToChange[i][j] -= matrixToChange[k][j] * me;
                }
            }
        }
        //возвращаем матрицу
        setMatrix(mCopy);
        inverseMatrix = matrixToChange;
        return inverseMatrix;
    }

    //печать матрицы
    public void printMatrix(double[][] matrix, int n) {
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                System.out.print(matrix[i][j] + " ");
            }
            System.out.println();
        }
    }

    //умножение матриц
    public double[][] multiplyMatrix(double[][] mA, double[][] mB, double[][] result, int n) {
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                for (int k = 0; k < n; k++) {
                    result[i][j] += mA[i][k] * mB[k][j];
                }
            }
        }
        return result;
    }

    //матричная норма
    public double matrNorm(double[][] a, int n) {

        double s;
        this.norm = 0.;

        for (int i = 0; i < n; i++) {
            s = 0;
            for (int j = 0; j < n; j++) {
                s += Math.abs(a[i][j]);
                if (s > norm) norm = s;
            }
        }

        return this.norm;
    }

    public void writeToFileHeader() throws IOException {
        try (PrintWriter pw = new PrintWriter(new FileOutputStream("t.txt", true))) {
            pw.printf("--------------------------------------------------------------------------------------------");
            pw.println();
            pw.print("|         a           |         b          |            R           |           Z          |");
            pw.println();
            pw.printf("--------------------------------------------------------------------------------------------");
            pw.println();
        }
    }

    public void writeToFile(double[][] m, double[][] mInv, double[][] aInv, double[][] r, int n, double alpha, double betta)
            throws IOException {
        try (PrintWriter out = new PrintWriter(new FileOutputStream("file.txt", true))) {

                double[] res = getRandZ(m, mInv, aInv, r, n);
                out.printf("| %14.10e |", alpha);
                out.printf(" %14.10e |", betta);
                out.printf(" %14.10e |", res[0]);
                out.printf(" %14.10e |", res[1]);
                out.printf(" %14.10e |", res[2]);
                out.printf("--------------------------------------------------------------------------------------------");
                out.println();


        }
    }

    public double[] getRandZ(double[][] m, double[][] mInv, double[][] aInv, double[][] r, int n) {
        double[] res = new double[3];
        for (int i = 0; i < n; i++) {
            r[i][i] -= 1;
        }
        multiplyMatrix(m, mInv, r, n);
        double norm = matrNorm(r, n);
        //System.out.print("|     " + norm + "      ");
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                r[i][j] = mInv[i][j] - aInv[i][j];
            }
        }
        res[0] = norm;
        norm = matrNorm(r, n);
        res[1] = norm;
        //System.out.print("|     " + norm + "      |");
        norm /= matrNorm(aInv, n);
        res[2] = norm;
        //System.out.print("     " + norm + "      |");
        return res;
    }


    public void getMaxElemAndCoeff() {
        double max = Double.MIN_VALUE;
        for (int i = 0; i < matrix.length; i++) {
            for (int j = 0; j < matrix[0].length; j++) {
                if (max < Math.abs(matrix[i][j])) {
                    maxLineElemNumber = i;
                    maxColumnElemNumber = j;
                    max = matrix[i][j];
                }
            }
        }
        maxElem = max;
    }

    public double[] getColumn(int j) {
        double[] column = new double[matrix[j].length];
        for (int i = 0; i < matrix[j].length; i++)
            column[i] = matrix[i][j];
        return column;
    }

    public double[] getLine(int i) {
        double[] line = new double[matrix[i].length];
        System.arraycopy(matrix[i], 0, line, 0, matrix[i].length);
        return line;
    }

    public void setLine(int k, double[] line, double[][] matrix) {
        for (int i = 0; i < line.length; i++) {
            matrix[k][i] = line[i];
        }
    }

    public double[][] getMatrix() {
        return matrix;
    }

    public void setMatrix(double[][] matrix) {
        for (int i = 0; i < getMatrix().length; i++) {
            for (int j = 0; j < getMatrix().length; j++) {
                this.matrix[i][j] = matrix[i][j];
            }
        }

    }

    public int getMaxLineElemNumber() {
        return maxLineElemNumber;
    }

    public void setMaxLineElemNumber(int maxLineElemNumber) {
        this.maxLineElemNumber = maxLineElemNumber;
    }

    public int getMaxColumnElemNumber() {
        return maxColumnElemNumber;
    }

    public void setMaxColumnElemNumber(int maxColumnElemNumber) {
        this.maxColumnElemNumber = maxColumnElemNumber;
    }

    public double getMaxElem() {
        double max = Double.MIN_VALUE;
        for (int i = 0; i < matrix.length; i++) {
            for (int j = 0; j < matrix[0].length; j++) {
                if (max < Math.abs(matrix[i][j])) {
                    maxElem = matrix[i][j];
                }
            }
        }
        return maxElem;
    }

    public void setMaxElem(double maxElem) {
        this.maxElem = maxElem;
    }


    @Override
    public String toString() {
        return "Matrix{" +
                "matrix=" + Arrays.toString(matrix) +
                '}';
    }
}

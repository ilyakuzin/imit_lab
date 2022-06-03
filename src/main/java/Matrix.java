import java.util.Arrays;

public class Matrix {

    private double[][] matrix;
    private int maxLineElemNumber;
    private int maxColumnElemNumber;
    private double maxElem;
    private double[][] inverseMatrix;

    public Matrix(double[][] matrix) {
        this.matrix = matrix;
    }

    //поиск обратной матрицы
    public double[][] getInverseMatrix() {
        double[][] matrixToChange = new double[matrix.length][matrix.length];
        double[][] mCopy = matrix;//сохраняем изначальную матрицу
        double me = maxElem;//коэффицент
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
        matrix = mCopy;
        inverseMatrix = matrixToChange;
        return inverseMatrix;
    }

    public void swapLines(int lineNumberToChange) {
        double[] line = getLine(lineNumberToChange);
        for (int j = 0; j < getMatrix().length; j++) {
            matrix[lineNumberToChange][j] = matrix[0][j];
            matrix[0][j] = line[j];
        }
    }

    public void swapColumn(int columnNumberToSwap) {
        double[] column = getColumn(columnNumberToSwap);
        for (int i = 0; i < getMatrix().length; i++) {
            matrix[i][columnNumberToSwap] = matrix[0][i];
            matrix[i][0] = column[i];
        }
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
        this.matrix = matrix;
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

    public void printMatrix(double[][] matrix) {
        for (int i = 0; i < matrix.length; i++) {
            for (int j = 0; j < matrix.length; j++) {
                System.out.print(matrix[i][j] + " ");
            }
            System.out.println();
        }
    }

    @Override
    public String toString() {
        return "Matrix{" +
                "matrix=" + Arrays.toString(matrix) +
                '}';
    }
}

import gen.Generator;

public class Main {

    public static void main(String[] args) {

        //размерность
        int n = 3;
        //матрицы
        double[][] m = new double[n][n];    //матрица
        double[][] mInv;                    //обратная матрица
        double[][] mI = new double[n][n];   //единичная

        //свои значения, некоторые значения матриц плохо считает, некоторые нормально, с это связано понять не могу
        //double[][] m = new double[][]{{3,5,6},{1,10,2},{1,3,4}};
        //double[][] m = new double[][]{{1,2,3},{44,33,22},{55,77,111}};
        //double[][] m = new double[][]{{1,0,0},{1,1,0},{0,0,2}}; //хорошо считает, простые матрицы

        Generator generator = new Generator();
        //нужное расскоментировать
        //generator.mygen(m, mInv, n, Math.pow(10, -10), 1, 1, 2,0, 1); // симметричная
        //generator.mygen(m, mInv, n, Math.pow(10, -14), 1, 1, 2, 1, 1); //простой структуры
        //generator.mygen(m, mInv, n, Math.pow(10, -14), 1, 0, 0, 2, 1); // жорданова клетка (n>2)

        //вывод матрицы
        System.out.println("|-------------------------------------------------------------------------------|");
        Matrix matrix = new Matrix(m);
        matrix.printMatrix(m, n);

        //вывод обратной матрицы
        System.out .println("|-------------------------------------------------------------------------------|");
        mInv = matrix.getInverseMatrix();
        matrix.printMatrix(mInv, n);

        //вывод умножения матриц
        System.out.println("|--------------------------------------------------------------------------------|");
        //generator.matr_mul(m, mInv, mI, n);
        //generator.print_matr(mI, n);
        matrix.printMatrix(matrix.multiplyMatrix(m, mInv, mI, n), n);


    }





}

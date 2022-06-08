import gen.Generator;

import java.io.IOException;

public class Main {

    public static void main(String[] args) throws IOException {

        //размерность
        int n = 10;
        //матрицы
        double[][] m = new double[n][n];       //матрица
        double[][] mInv = new double[n][n];     //обратная матрица
        double[][] aInv = new double[n][n];
        double[][] mI = new double[n][n];      //единичная
        double[][] r = new double[n][n];

        //свои значения, некоторые значения матриц плохо считает, некоторые нормально, с это связано понять не могу
        //double[][] m = new double[][]{{3,5,6},{1,10,2},{1,3,4}};
        //double[][] m = new double[][]{{1,2,3},{44,33,22},{55,77,111}};
        //double[][] m = new double[][]{{1,0,0},{1,1,0},{0,0,2}}; //хорошо считает, простые матрицы

        Generator generator = new Generator();
        //нужное расскоментировать
        //generator.mygen(m, mInv, n, Math.pow(10, -10), 1, 1, 2,0, 1); // симметричная
        generator.mygen(m, aInv, n, Math.pow(10, -14), 0.5, 1, 2, 1, 1); //простой структуры
        //generator.mygen(m, mInv, n, Math.pow(10, -14), 1, 0, 0, 2, 1); // жорданова клетка (n>2)
        //generator.genX()

        //вывод матрицы
        System.out.println("|-------------------------------------------------------------------------------|");
        Matrix matrix = new Matrix(m);
        //matrix.printMatrix(m, n);

        //вывод обратной матрицы
        System.out .println("|-------------------------------------------------------------------------------|");
        System.out.println("Обратная матрица");
        mInv = matrix.getInverseMatrix();
        matrix.printMatrix(mInv, n);

        //вывод обратной матрицы генератора
        System.out .println("|-------------------------------------------------------------------------------|");
        System.out.println("Обратная матрица генератора");
        matrix.printMatrix(aInv, n);

        //вывод умножения матриц
        System.out.println("|--------------------------------------------------------------------------------|");
        //generator.matr_mul(m, mInv, mI, n);
        //generator.print_matr(mI, n);
        matrix.printMatrix(matrix.multiplyMatrix(m, mInv, mI, n), n);

        //невязка и экперименты
        System.out.println("|--------------------------------------------------------------------------------|");
        System.out.println("Невязка и эксперименты");
        for (int i = 0; i < n; i++) {
            r[i][i] -= 1;
        }
        matrix.multiplyMatrix(m, mInv, r, n);
        double norm = matrix.matrNorm(r, n);
        System.out.println("|| R || = " + norm);
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                r[i][j] = mInv[i][j] - aInv[i][j];
            }
        }
        norm = matrix.matrNorm(r, n);
        System.out.println("||  Z  || = " + norm);
        norm /= matrix.matrNorm(aInv, n);
        System.out.println("||  Zeta   || = " + norm);

        //запись в файл
        //generator.writeHeadFile();
        //generator.writeFile();


        int i = 1;
        while (i <= 16) {
            double alpha = Math.pow(10, i);  //1;
            double betta = 1;               //Math.pow(10,i);
            generator.mygen(m, mInv, 10, alpha, betta, -1, 1, 2, 1);
            matrix.writeToFileHeader();
            matrix.writeToFile(m, mInv, aInv, r, n, alpha, betta);

            i++;
        }



    }





}

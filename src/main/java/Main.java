import gen.Generator;

public class Main {

    public static void main(String[] args) {

        int n = 10;
        double[][] m = new double[n][n];
        double [][] mInv = new double[n][n];
        Generator generator = new Generator();
        generator.mygen(m, mInv, n, Math.pow(10, -14), 1, 1, 2,0, 1); // симметричная
        //generator.mygen(m, mInv, n, Math.pow(10, -14), 1, 1, 2, 1, 1); //простой структуры
        //mygen(m, mInv, n, Math.pow(10, -14), 1, 0, 0, 2, 1); // жорданова клетка (n>2)
        System.out.println("|---------------------------------------------------------------------------------------------------------------------------------------------------------------------------|");
        //вывод матрицы
        Matrix matrix = new Matrix(m);
        double[][] invM = matrix.getInverseMatrix();
        matrix.printMatrix(invM);



    }





}

import java.util.Arrays;

public class Main {

    public static void main(String[] args) {

        double[][] m = {{1, 0, 0}, {0, 1, 0}, {1, 0, 1}};
        Matrix matrix = new Matrix(m);
        matrix.getMaxElemAndCoeff();
        Matrix iMatrix = new Matrix(matrix.getInverseMatrix());
        m = iMatrix.getMatrix();
        //вывод матрицы
        for (int i = 0; i < matrix.getMatrix().length; i++) {
            System.out.println(Arrays.toString(m[i]));
        }

    }





}

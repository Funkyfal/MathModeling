import java.text.DecimalFormat;
import java.util.Arrays;

import static java.lang.Math.*;

/**
 * вариант 2
 * 1) m = -3, s^2 = 16, sigma = 4
 */

public class Main {
    static class MCMGenerator {
        private final long beta;
        private final long M;
        private long state;

        public MCMGenerator(long beta, long seed, long M) {
            this.beta = beta;
            this.M = M;
            this.state = seed % M;
            if (this.state <= 0) {
                this.state += M;
            }
        }

        public MCMGenerator copy() {
            return new MCMGenerator(this.beta, this.state, this.M);
        }

        public long nextRaw() {
            state = ((beta * state) % M);
            return state;
        }

        public double nextDouble() {
            return ((double) nextRaw()) / ((double) M);
        }
    }
    
    static final DecimalFormat F = new DecimalFormat("0.00000");

    public static final int m1 = -3;
    public static final int sigma1 = 4;

    public static void main(String[] args) {
        final int n = 1000;
        final long M = 2147483648L; //модуль М

        long beta1 = 79507L;
        long a0_1 = 79507L;

        long beta2 = 65539L;
        long a0_2 = 256959681L;

        MCMGenerator firstGenerator = new MCMGenerator(beta1, a0_1, M);
        MCMGenerator secondGenerator = new MCMGenerator(beta2, a0_2, M);

        double[] normDistSample = generateNormDistribution(n, firstGenerator, secondGenerator);

        analyzeNormDist(normDistSample);
    }

    public static void analyzeNormDist(double[] normDistSample) {
        System.out.println("Нормальное одномерное распределение");
        System.out.println("    Среднее значение: " + F.format(mean(normDistSample)) + ", теоретически: " + m1);
        System.out.println("    Дисперсия: " + F.format(variance(normDistSample)) + ", теоретически: " + pow(sigma1, 2));

    }

    public static double[] generateNormDistribution(int n, MCMGenerator firstGenerator, MCMGenerator secondGenerator) {
        double[] dist = new double[n];
        for (int i = 0; i < n - 1; i += 2) {
            double a1 = firstGenerator.nextDouble();
            double a2 = secondGenerator.nextDouble();
            while (a1 == 0.0) {
                a1 = firstGenerator.nextDouble();
            }
            while (a2 == 0.0) {
                a2 = secondGenerator.nextDouble();
            }

            double n1 = sqrt(-2 * log(a1)) * cos(2 * PI * a2);
            double n2 = sqrt(-2 * log(a1)) * sin(2 * PI * a2);
            dist[i] = m1 + sigma1 * n1;
            dist[i + 1] = m1 + sigma1 * n2;
        }

        return dist;
    }

    public static double mean(double[] dist) {
        double sum = 0.0;
        for (double num : dist) {
            sum += num;
        }

        return sum / dist.length;
    }

    public static double variance(double[] dist) {
        double sum = 0.0;
        double mean = mean(dist);
        for (double num : dist) {
            sum += (num - mean) * (num - mean);
        }

        return sum / (dist.length - 1);
    }

//    public static boolean KolmogorovTest(double[] sample) {
//        double epsilon = 0.05;
//        double epsilonK = 1.358; // для epsilon=0.05
//        double delta = pow(epsilonK, -1) * (1 - epsilon);
//        double DnStatistics =
//    }
}


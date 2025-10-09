//вариант 2: a0=beta=79507, K=64, n=1000

import java.util.*;
import java.text.DecimalFormat;

public class Main {
    public static void main(String[] args) {
        final int n = 1000;
        final long M = 2147483648L; //модуль М

        //a0=beta=79507
        long beta1 = 79507L;
        long a0_1 = 79507L;

        //второй генератора для ММ
        long beta2 = 65539L;
        long a0_2 = 256959681L;

        int K = 64;

        MCMGenerator firstGenerator = new MCMGenerator(beta1, a0_1, M);
        MCMGenerator secondGenerator = new MCMGenerator(beta2, a0_2, M);

        //МКМ n образцов
        double[] sampleMCM = new double[n];
        for (int i = 0; i < n; i++) sampleMCM[i] = firstGenerator.nextDouble();

        //ММ n образцов
        MacLarenMarsaglia mm = new MacLarenMarsaglia(firstGenerator.copy(), secondGenerator, K);
        double[] sampleMM = new double[n];
        for (int i = 0; i < n; i++) sampleMM[i] = mm.nextDouble();

        DecimalFormat df = new DecimalFormat("#0.00000");

        System.out.println("результаты для n=" + n + "\n");

        System.out.println("1) МКМ (beta=" + beta1 + ", a0=" + a0_1 + ")");
        analyzeSample(sampleMCM, df);

        System.out.println("\n2) ММ (K=" + K + ", первый генератор - МКМ(beta=" + beta1 + ")," +
                " второй генератор - МКМ(beta=" + beta2 + "))");
        analyzeSample(sampleMM, df);
    }

    static void analyzeSample(double[] sample, DecimalFormat df) {
        int n = sample.length;
        double mean = mean(sample);
        double var = variance(sample, mean);
        System.out.println("среднее значение (матожидание): " + df.format(mean) + "  (теор. 0.5)");
        System.out.println("дисперсия: " + df.format(var) + "  (теор. ~0.083333)");

        //критерий колмогорова
        KSTest.Result ks = KSTest.kolmogorovSmirnovTest(sample);
        System.out.println("\nКолмогоров-Смирнов:");
        double alpha = 0.05;
        double alphaK = 1.358; // для alpha=0.05
        double criticalD = alphaK / Math.sqrt(n);
        System.out.println(" критическое D (α=" + alpha + ") = " + df.format(criticalD));
        System.out.println(" статистика: " + ks.D);
        System.out.println(" решение по критическому значению: " + (ks.D > criticalD ? "отвергаем H0" : "не отвергаем H0"));

        // χ^2 пирсона разбиение на 10 равных ячеек
        int bins = 10;
        ChiSquareTest.Result ch = ChiSquareTest.chiSquareTestUniform(sample, bins);
        int dfChi = (bins - 1); //степень свободы
        System.out.println("\nХ квадрат:");
        double chi2Crit = 16.9190; // табличное значение для df=9, alpha=0.05
        System.out.println(" Х квадрат критическое (df=" + dfChi + ", α=" + alpha + ") = " + df.format(chi2Crit));
        System.out.println(" Х квадрат статистика: " + ch.chi2);
        System.out.println(" Решение по критическому значению: " + (ch.chi2 > chi2Crit ? "отвергаем H0" : "не отвергаем H0"));
    }

    static double mean(double[] a) {
        double s = 0;
        for (double v : a) {
            s += v;
        }
        return s / a.length;
    }

    static double variance(double[] a, double mean) {
        double s = 0;
        for (double v : a) {
            s += (v - mean) * (v - mean);
        }
        return s / (a.length - 1);
    }

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

    static class MacLarenMarsaglia {
        private final MCMGenerator firstGenerator;
        private final MCMGenerator secondGenerator;
        private final double[] table;
        private final int K;

        public MacLarenMarsaglia(MCMGenerator firstGenerator, MCMGenerator secondGenerator, int K) {
            this.firstGenerator = firstGenerator;
            this.secondGenerator = secondGenerator;
            this.K = K;
            this.table = new double[K];
            for (int i = 0; i < K; i++) table[i] = firstGenerator.nextDouble(); //сначала надо заполнить таблицу
        }

        public double nextDouble() {
            double cb = secondGenerator.nextDouble();//нужно для выбора индекса
            int idx = (int) (Math.floor(cb * K));
            if (idx < 0) {
                idx = 0;
            }
            if (idx >= K) {
                idx = K - 1;
            }
            double result = table[idx];
            table[idx] = firstGenerator.nextDouble(); //заносим в таблицу новое значение
            return result;
        }
    }

    //колмогоров-смирнов для того, чтобы проверить на равномерность
    static class KSTest {
        static class Result {
            double D;
            double sqrtNTimesD;
        }

        static Result kolmogorovSmirnovTest(double[] sample) {
            int n = sample.length;
            double[] s = Arrays.copyOf(sample, n);
            Arrays.sort(s);
            double D = 0.0;
            for (int i = 0; i < n; i++) {
                double Fi = (i + 1) / (double) n;
                double Fim1 = (i) / (double) n;
                double xi = s[i];
                D = Math.max(D, Math.max(Math.abs(Fi - xi), Math.abs(Fim1 - xi)));
            }

            double y = Math.sqrt(n) * D;
            Result r = new Result();
            r.D = D;
            r.sqrtNTimesD = y;
            return r;
        }
    }

    //Хи в квадрате
    static class ChiSquareTest {
        static class Result {
            double chi2;
            int df;
        }

        //разбиваем на бины
        static Result chiSquareTestUniform(double[] sample, int bins) {
            int n = sample.length;
            int[] obs = new int[bins];
            for (double v : sample) {
                int idx = (int) Math.floor(v * bins);
                if (idx < 0) {
                    idx = 0;
                }
                if (idx >= bins) {
                    idx = bins - 1;
                }
                obs[idx]++;
            }
            double expected = n / (double) bins;
            double chi2 = 0.0;
            for (int i = 0; i < bins; i++) {
                double diff = obs[i] - expected;
                chi2 += diff * diff / expected;
            }
            int df = bins - 1;
            Result r = new Result();
            r.chi2 = chi2;
            r.df = df;
            return r;
        }
    }
}

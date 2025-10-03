// Lab1_BSV.java
// Лабораторная работа №1 — моделирование БСВ
// Вариант 2: a0 = beta = 79507, K = 64, n = 1000
import java.util.*;
import java.text.DecimalFormat;

public class Main {
    public static void main(String[] args) {
        final int n = 1000;
        final long M = 2147483648L; // 2^31 - 1 (типичный модуль для МКМ в учебнике). См. пособие. :contentReference[oaicite:3]{index=3}

        // Вариант 2: a0 = beta = 79507, K=64
        long beta1 = 79507L;
        long a0_1  = 79507L;

        // Второй генератор для Макларена-Марсальи возьмём образец из пособия (пример D2). См. пособие. :contentReference[oaicite:4]{index=4}
        long beta2 = 65539L;
        long a0_2  = 256959681L;

        int K = 64;

        MCMGenerator g1 = new MCMGenerator(beta1, a0_1, M);
        MCMGenerator g2 = new MCMGenerator(beta2, a0_2, M);

        // 1) Сгенерировать n реализаций МКМ-датчика (g1)
        double[] sampleMCM = new double[n];
        for (int i = 0; i < n; i++) sampleMCM[i] = g1.nextDouble();

        // 2) Сгенерировать n реализаций методом Макларена-Марсальи:
        MacLarenMarsaglia mm = new MacLarenMarsaglia(g1.copy(), g2, K); // g1.copy() — один из простейших датчиков (МКМ из п.1)
        double[] sampleMM = new double[n];
        for (int i = 0; i < n; i++) sampleMM[i] = mm.nextDouble();

        DecimalFormat df = new DecimalFormat("#0.00000");

        System.out.println("=== Результаты моделирования (n=" + n + ") ===\n");

        System.out.println("== 1) МКМ-датчик (beta=" + beta1 + ", a0=" + a0_1 + ") ==");
        analyzeSample(sampleMCM, df);

        System.out.println("\n== 2) Макларен-Марсальи (K=" + K + ", D1=МКМ(beta=" + beta1 + "), D2=МКМ(beta=" + beta2 + ")) ==");
        analyzeSample(sampleMM, df);
    }

    static void analyzeSample(double[] sample, DecimalFormat df) {
        int n = sample.length;
        double mean = mean(sample);
        double var = variance(sample, mean);
        System.out.println("Среднее (эст.): " + df.format(mean) + "  (теоретически 0.5)");
        System.out.println("Дисперсия (эст.): " + df.format(var) + "  (теоретически ~0.083333)");

        // Критерий Колмогорова
        KSTest.Result ks = KSTest.kolmogorovSmirnovTest(sample);
        System.out.println("\nKolmogorov-Smirnov:");
        double alpha = 0.05;
        double K_alpha = 1.358; // для alpha=0.05 (двусторонний тест)
        double Dcrit = K_alpha / Math.sqrt(n);
        System.out.println(" K-S критическое D (α=" + alpha + ") = " + df.format(Dcrit));
        System.out.println(" K-S статистика: " + ks.D);
        System.out.println(" Решение по критическому значению: " + (ks.D > Dcrit ? "отвергаем H0" : "не отвергаем H0"));

        // χ^2 Пирсона — разбиение на 10 равных ячеек (ожидание = n/10)
        int bins = 10;
        ChiSquareTest.Result ch = ChiSquareTest.chiSquareTestUniform(sample, bins);
        int dfChi = (bins - 1);
        System.out.println("\nChi-square:");
        double chi2Crit = 16.9190; // табличное значение для df=9, alpha=0.05
        System.out.println(" Chi-square критическое (df=" + dfChi + ", α=" + alpha + ") = " + df.format(chi2Crit));
        System.out.println(" X² статистика: " + ch.chi2);
        System.out.println(" Решение по критическому значению: " + (ch.chi2 > chi2Crit ? "отвергаем H0" : "не отвергаем H0"));
    }

    static double mean(double[] a) {
        double s = 0;
        for (double v : a) s += v;
        return s / a.length;
    }
    static double variance(double[] a, double mean) {
        double s = 0;
        for (double v : a) s += (v - mean)*(v - mean);
        return s / (a.length - 1);
    }

    // ---------- MCM generator ----------
    static class MCMGenerator {
        private final long beta;
        private final long M;
        private long state;

        public MCMGenerator(long beta, long seed, long M) {
            this.beta = beta;
            this.M = M;
            this.state = seed % M;
            if (this.state <= 0) this.state += M; // стартовое значение должно быть в 1..M-1
        }

        public MCMGenerator copy() {
            return new MCMGenerator(this.beta, this.state, this.M);
        }
        public long nextRaw() {
            state = ( (beta * state) % M );
            return state;
        }
        public double nextDouble() {
            return ((double) nextRaw()) / ((double) M);
        }
    }

    // ---------- MacLaren-Marsaglia ----------
    static class MacLarenMarsaglia {
        private final MCMGenerator A;
        private final MCMGenerator B;
        private final double[] V;
        private final int K;
        public MacLarenMarsaglia(MCMGenerator A, MCMGenerator B, int K) {
            // A и B должны быть независимы (в нашем случае используются разные betas/seeds)
            this.A = A;
            this.B = B;
            this.K = K;
            this.V = new double[K];
            for (int i = 0; i < K; i++) V[i] = A.nextDouble(); // первоначальное заполнение таблицы
        }
        public double nextDouble() {
            double cb = B.nextDouble(); // для выбора индекса
            int idx = (int) (Math.floor(cb * K));
            if (idx < 0) idx = 0;
            if (idx >= K) idx = K - 1;
            double result = V[idx];
            V[idx] = A.nextDouble(); // обновление таблицы
            return result;
        }
    }

    // ---------- Kolmogorov-Smirnov ----------
    static class KSTest {
        static class Result { double D; double sqrtNtimesD; }
        // Для проверки равномерности на [0,1)
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
            double p = ksPvalue(y);
            Result r = new Result();
            r.D = D;
            r.sqrtNtimesD = y;
            return r;
        }

        // Kolmogorov distribution tail approximation:
        // Q(y) = 2 * sum_{k=1..inf} (-1)^{k-1} exp(-2 k^2 y^2)
        static double ksPvalue(double y) {
            if (y <= 0) return 1.0;
            double sum = 0.0;
            double term;
            int k = 1;
            do {
                term = 2.0 * Math.pow(-1.0, k-1) * Math.exp(-2.0 * k * k * y * y);
                sum += term;
                k++;
            } while (Math.abs(term) > 1e-12 && k < 200);
            double p = sum;
            if (p < 0) p = 0.0;
            if (p > 1) p = 1.0;
            return p;
        }
    }

    // ---------- Chi-square test for uniformity ----------
    static class ChiSquareTest {
        static class Result { double chi2; int df; }
        // Разбиваем [0,1) на 'bins' равных интервалов
        static Result chiSquareTestUniform(double[] sample, int bins) {
            int n = sample.length;
            int[] obs = new int[bins];
            for (double v : sample) {
                int idx = (int) Math.floor(v * bins);
                if (idx < 0) idx = 0;
                if (idx >= bins) idx = bins - 1;
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

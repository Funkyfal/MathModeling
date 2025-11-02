import java.text.DecimalFormat;
import java.util.Arrays;

/**
 * Вариант 2
 * <p>1) m = -3, s^2 = 16, sigma = 4</p>
 * <p>2) m = 0, s^2 = 4</p>
 * <p>3) a = 1, b = 2</p>
 */
public class Main {
    static final DecimalFormat F = new DecimalFormat("0.00000");

    public static final double mNorm = -3.0;
    public static final double s2Norm = 16.0;
    public static final double sigmaNorm = Math.sqrt(s2Norm);

    public static final double mLog = 0.0;
    public static final double s2Log = 4.0;
    public static final double sigmaLog = Math.sqrt(s2Log);

    public static final double aCauchy = 1.0;
    public static final double bCauchy = 2.0;

    public static void main(String[] args) {
        final int n = 1000;

        final long M = 2147483648L;
        long beta1 = 79507L;
        long a0_1 = 79507L;
        long beta2 = 65539L;
        long a0_2 = 256959681L;

        MCMGenerator gen1 = new MCMGenerator(beta1, a0_1, M);
        MCMGenerator gen2 = new MCMGenerator(beta2, a0_2, M);

        double[] sampleNorm = generateNormal(n, gen1, gen2, mNorm, sigmaNorm);
        analyzeSample("Нормальное N(-3,16)", sampleNorm, "norm", new double[]{mNorm, s2Norm});

        MCMGenerator g1ln = new MCMGenerator(beta1, a0_1, M);
        MCMGenerator g2ln = new MCMGenerator(beta2, a0_2, M);
        double[] sampleLogN = generateLogNormal(n, g1ln, g2ln, mLog, s2Log);
        analyzeSample("Логнормальное LN(0,4)", sampleLogN, "lognorm", new double[]{mLog, s2Log});

        MCMGenerator gc = new MCMGenerator(beta1, a0_1, M);
        double[] sampleCauchy = generateCauchy(n, gc, aCauchy, bCauchy);
        analyzeSample("Коши C(1,2)", sampleCauchy, "cauchy", new double[]{aCauchy, bCauchy});

        final int REPS = 2000;
        System.out.println("\n--- Монте-Карло: оценка вероятности ошибки I рода (повторения = " + REPS + ") ---");

        double type1_ks_norm = estimateTypeI_KS("norm", new double[]{mNorm, s2Norm}, n, REPS, beta1 + 10, a0_1 + 10, M);
        System.out.println("Колмогоров ошибка первого рода (норм)  ≈ " + F.format(type1_ks_norm));

        double type1_chi_norm = estimateTypeI_Chi2("norm", new double[]{mNorm, s2Norm}, n, REPS, beta1 + 20, a0_1 + 20, M);
        System.out.println("Х²ошибка первого рода (норм) ≈ " + F.format(type1_chi_norm));

        double type1_ks_logn = estimateTypeI_KS("lognorm", new double[]{mLog, s2Log}, n, REPS, beta1 + 30, a0_1 + 30, M);
        System.out.println("Колмогоров ошибка первого рода (логнорм)  ≈ " + F.format(type1_ks_logn));

        double type1_chi_logn = estimateTypeI_Chi2("lognorm", new double[]{mLog, s2Log}, n, REPS, beta1 + 40, a0_1 + 40, M);
        System.out.println("Х² ошибка первого рода (логнорм) ≈ " + F.format(type1_chi_logn));

        double type1_ks_cauchy = estimateTypeI_KS("cauchy", new double[]{aCauchy, bCauchy}, n, REPS, beta1 + 50, a0_1 + 50, M);
        System.out.println("Колмогоров ошибка первого рода (коши)   ≈ " + F.format(type1_ks_cauchy));

        double type1_chi_cauchy = estimateTypeI_Chi2("cauchy", new double[]{aCauchy, bCauchy}, n, REPS, beta1 + 60, a0_1 + 60, M);
        System.out.println("Х² ошибка первого рода (коши) ≈ " + F.format(type1_chi_cauchy));
    }

    public static void analyzeSample(String title, double[] sample, String distName, double[] params) {
        System.out.println("\n=== " + title + " ===");
        System.out.println("  n = " + sample.length);
        System.out.println("  матожидание: " + F.format(mean(sample)));
        System.out.println("  дисперсия: " + F.format(variance(sample)));

        KSTest.Result ksres = KSTest.kolmogorovSmirnovTest(sample, distName, params);
        System.out.println("  Колмогоров: Dn = " + F.format(ksres.D) + ", sqrt(n)*Dn = " + F.format(ksres.sqrtNTimesD)
                + ", принимаем гипотезу H0? " + ksres.acceptH0);

        // Chi2
        Chi2Test.Result chires = Chi2Test.chi2Test(sample, distName, params, 10);
        System.out.println("  Х²: статистика = " + F.format(chires.chi2) + ", СС = " + chires.df
                + ", критическое(0.95) = " + F.format(chires.chi2Crit) + ", принимаем гипотезу H0? " + chires.acceptH0);
    }

    public static double[] generateNormal(int n, MCMGenerator g1, MCMGenerator g2, double m, double sigma) {
        double[] sample = new double[n];
        for (int i = 0; i < n; i += 2) {
            double u1 = g1.nextDouble();
            while (u1 <= 0.0)
                u1 = g1.nextDouble();
            double u2 = g2.nextDouble();
            while (u2 <= 0.0)
                u2 = g2.nextDouble();

            double r = Math.sqrt(-2.0 * Math.log(u1));
            double z1 = r * Math.cos(2.0 * Math.PI * u2);
            double z2 = r * Math.sin(2.0 * Math.PI * u2);

            sample[i] = m + sigma * z1;
            if (i + 1 < n) sample[i + 1] = m + sigma * z2;
        }
        return sample;
    }

    public static double[] generateLogNormal(int n, MCMGenerator g1, MCMGenerator g2, double mLog, double s2Log) {
        double[] z = generateNormal(n, g1, g2, mLog, Math.sqrt(s2Log));
        double[] sample = new double[n];
        for (int i = 0; i < n; i++)
            sample[i] = Math.exp(z[i]);
        return sample;
    }

    public static double[] generateCauchy(int n, MCMGenerator g, double a, double b) {
        double[] sample = new double[n];
        for (int i = 0; i < n; i++) {
            double u = g.nextDouble();
            while (u <= 0.0 || u >= 1.0)
                u = g.nextDouble();
            sample[i] = a + b * Math.tan(Math.PI * (u - 0.5));
        }
        return sample;
    }

    public static double mean(double[] dist) {
        double sum = 0.0;
        for (double v : dist) sum += v;
        return sum / dist.length;
    }

    public static double variance(double[] dist) {
        double mean = mean(dist);
        double sum = 0.0;
        for (double v : dist) sum += (v - mean) * (v - mean);
        return sum / (dist.length - 1);
    }

    static class KSTest {
        static class Result {
            double D;
            double sqrtNTimesD;
            boolean acceptH0;
        }

        static Result kolmogorovSmirnovTest(double[] sample, String distName, double[] params) {
            int n = sample.length;
            double[] u = new double[n];
            for (int i = 0; i < n; i++) {
                u[i] = theoreticalCdf(sample[i], distName, params);
            }
            Arrays.sort(u);

            double Dplus = 0.0, Dminus = 0.0;
            for (int i = 1; i <= n; i++) {
                double Ui = u[i - 1];
                double EmpUp = (double) i / n;
                double EmpLow = (double) (i - 1) / n;
                Dplus = Math.max(Dplus, EmpUp - Ui);
                Dminus = Math.max(Dminus, Ui - EmpLow);
            }
            double D = Math.max(Dplus, Dminus);
            double sqrtNTimesD = Math.sqrt(n) * D;

            double Kalpha = 1.358; // таблица для alpha=0.05
            double Dcrit = Kalpha / Math.sqrt(n);
            boolean accept = D < Dcrit;

            Result r = new Result();
            r.D = D;
            r.sqrtNTimesD = sqrtNTimesD;
            r.acceptH0 = accept;
            return r;
        }
    }

    static class Chi2Test {
        static class Result {
            double chi2;
            int df;
            double chi2Crit;
            boolean acceptH0;
            int[] obs;
            double[] expected;
        }

        static Result chi2Test(double[] sample, String distName, double[] params, int k) {
            int n = sample.length;
            double[] probs = new double[k + 1];
            for (int i = 0; i <= k; i++)
                probs[i] = (double) i / k;

            double[] edges = new double[k + 1];
            for (int i = 0; i <= k; i++) {
                double p = probs[i];
                if (p <= 0.0) edges[i] = -Double.MAX_VALUE;
                else if (p >= 1.0) edges[i] = Double.MAX_VALUE;
                else edges[i] = inverseCdf(p, distName, params);
            }

            int[] obs = new int[k];
            for (double x : sample) {
                int j = Arrays.binarySearch(edges, x);
                if (j >= 0) {
                    j = Math.min(j, k - 1);
                } else {
                    j = -j - 2;
                    if (j < 0) j = 0;
                    if (j >= k) j = k - 1;
                }
                obs[j]++;
            }

            double expectedCount = (double) n / k;
            double chi2 = 0.0;
            for (int i = 0; i < k; i++) {
                double diff = obs[i] - expectedCount;
                chi2 += diff * diff / expectedCount;
            }

            int df = k - 1;
            double chi2Crit = 16.919; // для k=10
            boolean accept = chi2 < chi2Crit;

            Result r = new Result();
            r.chi2 = chi2;
            r.df = df;
            r.chi2Crit = chi2Crit;
            r.acceptH0 = accept;
            r.obs = obs;
            r.expected = new double[k];
            Arrays.fill(r.expected, expectedCount);
            return r;
        }
    }

    private static double theoreticalCdf(double x, String distName, double[] params) {
        switch (distName) {
            case "norm": {
                double m = params[0];
                double s2 = params[1];
                double sigma = Math.sqrt(s2);
                return normalCdf(x, m, sigma);
            }
            case "lognorm": {
                double m = params[0];
                double s2 = params[1];
                double sigma = Math.sqrt(s2);
                if (x <= 0.0) return 0.0;
                return normalCdf(Math.log(x), m, sigma);
            }
            case "cauchy": {
                double a = params[0], b = params[1];
                return 0.5 + (1.0 / Math.PI) * Math.atan((x - a) / b);
            }
            default:
                throw new IllegalArgumentException("Unknown distName: " + distName);
        }
    }

    private static double inverseCdf(double prob, String distName, double[] params) {
        if ("cauchy".equals(distName)) {
            double a = params[0], b = params[1];
            return a + b * Math.tan(Math.PI * (prob - 0.5));
        } else if ("norm".equals(distName)) {
            double m = params[0];
            double s2 = params[1];
            double sigma = Math.sqrt(s2);
            return inverseNormalCdf(prob, m, sigma);
        } else if ("lognorm".equals(distName)) {
            double m = params[0];
            double s2 = params[1];
            double sigma = Math.sqrt(s2);
            double q = inverseNormalCdf(prob, m, sigma);
            return Math.exp(q);
        } else {
            throw new IllegalArgumentException("Unknown distName: " + distName);
        }
    }

    private static double inverseNormalCdf(double p, double mean, double sigma) {
        if (p <= 0.0) {
            return -Double.MAX_VALUE;
        }
        if (p >= 1.0) {
            return Double.MAX_VALUE;
        }
        double low = mean - 10 * sigma;
        double high = mean + 10 * sigma;
        for (int iter = 0; iter < 100; iter++) {
            double mid = 0.5 * (low + high);
            double c = normalCdf(mid, mean, sigma);
            if (c > p) {
                high = mid;
            } else {
                low = mid;
            }
        }
        return 0.5 * (low + high);
    }

    private static double normalCdf(double x, double mean, double sigma) {
        double z = (x - mean) / (sigma * Math.sqrt(2.0));
        return 0.5 * (1.0 + erf(z));
    }

    // аппроксимация erf
    private static double erf(double z) {
        double t = 1.0 / (1.0 + 0.5 * Math.abs(z));
        double tau = t * Math.exp(-z * z - 1.26551223 +
                t * (1.00002368 +
                        t * (0.37409196 +
                                t * (0.09678418 +
                                        t * (-0.18628806 +
                                                t * (0.27886807 +
                                                        t * (-1.13520398 +
                                                                t * (1.48851587 +
                                                                        t * (-0.82215223 +
                                                                                t * 0.17087277)))))))));
        return z >= 0 ? 1.0 - tau : tau - 1.0;
    }

    private static double estimateTypeI_KS(String distName, double[] params, int n, int reps, long seedBase1, long seedBase2, long M) {
        int rejections = 0;
        for (int r = 0; r < reps; r++) {
            MCMGenerator g1 = new MCMGenerator(79507L, seedBase1 + r, M);
            MCMGenerator g2 = new MCMGenerator(65539L, seedBase2 + r, M);

            double[] sample;
            if ("norm".equals(distName)) {
                sample = generateNormal(n, g1, g2, params[0], Math.sqrt(params[1]));
            } else if ("lognorm".equals(distName)) {
                sample = generateLogNormal(n, g1, g2, params[0], params[1]);
            } else if ("cauchy".equals(distName)) {
                sample = generateCauchy(n, g1, params[0], params[1]);
            } else throw new IllegalArgumentException();

            KSTest.Result res = KSTest.kolmogorovSmirnovTest(sample, distName, params);
            if (!res.acceptH0) rejections++;
        }
        return (double) rejections / reps;
    }

    private static double estimateTypeI_Chi2(String distName, double[] params, int n, int reps, long seedBase1, long seedBase2, long M) {
        int rejections = 0;
        for (int r = 0; r < reps; r++) {
            MCMGenerator g1 = new MCMGenerator(79507L, seedBase1 + r, M);
            MCMGenerator g2 = new MCMGenerator(65539L, seedBase2 + r, M);

            double[] sample;
            if ("norm".equals(distName)) {
                sample = generateNormal(n, g1, g2, params[0], Math.sqrt(params[1]));
            } else if ("lognorm".equals(distName)) {
                sample = generateLogNormal(n, g1, g2, params[0], params[1]);
            } else if ("cauchy".equals(distName)) {
                sample = generateCauchy(n, g1, params[0], params[1]);
            } else throw new IllegalArgumentException();

            Chi2Test.Result cres = Chi2Test.chi2Test(sample, distName, params, 10);
            if (!cres.acceptH0) rejections++;
        }
        return (double) rejections / reps;
    }
}

class MCMGenerator {
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

    public long nextRaw() {
        state = ((beta * state) % M);
        if (state <= 0)
            state += M;
        return state;
    }

    public double nextDouble() {
        return ((double) nextRaw()) / ((double) M);
    }
}

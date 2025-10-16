// монте-карло: java Lab2DiscreteNoPValue --mc 500

import java.util.*;
import java.text.DecimalFormat;

public class Main {
    static final DecimalFormat F = new DecimalFormat("0.00000");

    // критические значения Хи квадрат для альфа = 0.05, df = 1..50
    static final double[] CHI2_CRIT_95 = {
            0.0,
            3.841458820694124,
            5.991464547107979,
            7.814727903251179,
            9.487729036781154,
            11.070497693516351,
            12.591587243743977,
            14.067140449340169,
            15.50731305586545,
            16.918977604620448,
            18.307038053275146,
            19.67513657172878,
            21.02606981748379,
            22.36203208523671,
            23.68479109926933,
            24.99579013972857,
            26.29622760432735,
            27.587111572777997,
            28.86929910618511,
            30.143527233257095,
            31.410432844059243,
            32.67056460143622,
            33.92445196154039,
            35.17246397113603,
            36.4150293853458,
            37.65246635017661,
            38.88501775351129,
            40.11299851614626,
            41.33648085985938,
            42.55589263286622,
            43.77300142962844,
            44.98529962098502,
            46.19409214294103,
            47.39961439341192,
            48.60203406748361,
            49.80151291732866,
            50.99821960095802,
            52.19231101606287,
            53.38393672238559,
            54.57324603920106,
            55.76037415395297,
            56.945451017375,
            58.128589006431,
            59.309893186018,
            60.489462225262,
            61.667387742617,
            62.843754294396,
            64.018641422735,
            65.192123026177,
            66.364267137463,
            67.535133070786
    };

    public static void main(String[] args) {
        int n = 1000;
        double alpha = 0.05;

        double pBern = 0.5;
        int rNegBin = 5;
        double pNegBin = 0.25;

        if (args.length >= 2 && args[0].equals("--mc")) {
            int Nexp = Integer.parseInt(args[1]);
            monteCarloEstimateTypeI(n, pBern, rNegBin, pNegBin, alpha, Nexp);
            return;
        }

        Random rnd = new Random();

        int[] bernSample = generateBernoulliSample(n, pBern, rnd);
        int[] nbSample = generateNegBinomialSample(n, rNegBin, pNegBin, rnd);

        System.out.println("однократный прогон при n = " + n + "\n");

        System.out.println("Бернулли (1, p=" + pBern + "):");
        analyzeAndTestBernoulli(bernSample, pBern, alpha);

        System.out.println("\nОтрицательное биномиальное (r=" + rNegBin + ", p=" + pNegBin + "):");
        analyzeAndTestNegBin(nbSample, rNegBin, pNegBin, alpha);
    }

    //генерация
    static int[] generateBernoulliSample(int n, double p, Random rng) {
        int[] s = new int[n];
        for (int i = 0; i < n; i++) {
            s[i] = (rng.nextDouble() < p) ? 1 : 0;
        }
        return s;
    }

    static int[] generateNegBinomialSample(int n, int r, double p, Random rng) {
        int[] s = new int[n];
        for (int i = 0; i < n; i++) {
            int successes = 0;
            int failures = 0;
            while (successes < r) {
                if (rng.nextDouble() < p) {
                    successes++;
                } else {
                    failures++;
                }
            }
            s[i] = failures;
        }
        return s;
    }

    //статистика
    static double mean(int[] x) {
        double sum = 0;
        for (int v : x) {
            sum += v;
        }
        return sum / x.length;
    }

    static double sampleVariance(int[] x) {
        int n = x.length;
        double m = mean(x);
        double s = 0;
        for (int v : x) {
            s += (v - m) * (v - m);
        }
        return s / (n - 1);
    }

    //бернулли
    static void analyzeAndTestBernoulli(int[] sample, double p, double alpha) {
        int n = sample.length;
        double m = mean(sample);
        double var = sampleVariance(sample);
        double theorMean = p;
        double theorVar = p * (1 - p);
        System.out.println("среднее значение (матожидание) = " + F.format(m) + " (теоретически " + F.format(theorMean) + ")");
        System.out.println("дисперсия  = " + F.format(var) + " (теоретически " + F.format(theorVar) + ")");

        int c0 = 0;
        for (int v : sample) {
            if (v == 0) c0++;
        }
        int c1 = n - c0;
        double exp0 = n * (1 - p);
        double exp1 = n * p;
        double chi2 = (c0 - exp0) * (c0 - exp0) / exp0 + (c1 - exp1) * (c1 - exp1) / exp1;

        int df = 1; //2 бина => степень свободы = 1
        double crit = chi2Critical(df, alpha);
        System.out.println("\nХи квадрат (бернулли): статистика = " + F.format(chi2) + ", критическое(" + alpha + ", СС=" + df + ") = " + F.format(crit));
        System.out.println("решение: " + (chi2 <= crit ? "принимаем гипотезу (подходит)" : "отклоняем гипотезу"));
    }

    //негативное биномиальное
    static void analyzeAndTestNegBin(int[] sample, int r, double p, double alpha) {
        int n = sample.length;
        double m = mean(sample);
        double var = sampleVariance(sample);
        double theorMean = r * (1 - p) / p;
        double theorVar = r * (1 - p) / (p * p);
        System.out.println("среднее значение (матожидание) = " + F.format(m) + " (теоретически " + F.format(theorMean) + ")");
        System.out.println("дисперсия  = " + F.format(var) + " (теоретически " + F.format(theorVar) + ")");

        List<Double> pmf = new ArrayList<>();
        double p0 = Math.pow(p, r);
        pmf.add(p0);
        double tail = p0;
        int i = 0;
        int maxIter = 10000;
        while (tail < 1.0 - 1e-12 && i < maxIter) {
            double prev = pmf.get(i);
            double next = prev * ((double) (i + r) / (i + 1)) * (1 - p);
            pmf.add(next);
            tail += next;
            i++;
            if (i > 5000) {
                break;
            }
        }

        double sumProb = 0;
        for (double q : pmf) {
            sumProb += q;
        }
        if (sumProb < 1.0) {
            double last = pmf.get(pmf.size() - 1);
            pmf.set(pmf.size() - 1, last + (1.0 - sumProb));
        }

        List<Bin> bins = makeBinsFromPmf(pmf, n, 5.0);

        int[] obs = new int[bins.size()];
        int maxIndex = pmf.size() - 1;
        int[] valToBin = new int[maxIndex + 1];
        Arrays.fill(valToBin, -1);
        for (int bi = 0; bi < bins.size(); bi++) {
            for (int v = bins.get(bi).from; v <= bins.get(bi).to; v++) {
                if (v <= maxIndex) {
                    valToBin[v] = bi;
                }
            }
        }
        for (int v : sample) {
            int bi;
            if (v <= maxIndex) {
                bi = valToBin[v];
                if (bi == -1) {
                    bi = bins.size() - 1;
                }
            } else {
                bi = bins.size() - 1;
            }
            obs[bi]++;
        }

        double chi2 = 0;
        for (int bi = 0; bi < bins.size(); bi++) {
            double expected = bins.get(bi).prob * n;
            double o = obs[bi];
            double diff = o - expected;
            chi2 += diff * diff / expected;
        }
        int df = Math.max(1, bins.size() - 1);//степени свободы
        double crit = chi2Critical(df, alpha);

        System.out.println("\nХи квадрат (отрицательный биномиальный): статистика = " + F.format(chi2) + ", СС = " + df + ", критическое(" + alpha + ") = " + F.format(crit));
        System.out.println("бинов использовано (от..до) и предполагаемое количество:");
        for (int bi = 0; bi < bins.size(); bi++) {
            System.out.println("  [" + bins.get(bi).from + ".." + bins.get(bi).to + "]  ожидаем=" + F.format(bins.get(bi).prob * n) + "  наблюдаем=" + obs[bi]);
        }
        System.out.println("Решение: " + (chi2 <= crit ? "принимаем гипотезу (подходит)" : "отклоняем гипотезу"));
    }

    //помогает сделать наполнение бинов >= минимального количества
    static class Bin {
        int from, to;
        double prob;

        Bin(int f, int t, double p) {
            from = f;
            to = t;
            prob = p;
        }
    }

    static List<Bin> makeBinsFromPmf(List<Double> pmf, int n, double minExpected) {
        List<Bin> bins = new ArrayList<>();
        int i = 0;
        int M = pmf.size();
        while (i < M) {
            double accProb = 0;
            int start = i;
            while (i < M && accProb * n < minExpected) {
                accProb += pmf.get(i);
                i++;
            }
            if (i >= M) {
                bins.add(new Bin(start, M - 1, accProb));
                break;
            } else {
                bins.add(new Bin(start, i - 1, accProb));
            }
        }
        if (bins.size() >= 2) {
            Bin last = bins.get(bins.size() - 1);
            if (last.prob * n < minExpected) {
                Bin prev = bins.get(bins.size() - 2);
                Bin merged = new Bin(prev.from, last.to, prev.prob + last.prob);
                bins.set(bins.size() - 2, merged);
                bins.remove(bins.size() - 1);
            }
        }
        return bins;
    }

    //монте-карло
    static void monteCarloEstimateTypeI(int n, double pBern, int rNegBin, double pNegBin, double alpha, int Nexp) {
        Random rng = new Random();
        int rejectBern = 0;
        int rejectNegBin = 0;
        for (int t = 0; t < Nexp; t++) {
            int[] b = generateBernoulliSample(n, pBern, rng);
            int[] nb = generateNegBinomialSample(n, rNegBin, pNegBin, rng);
            if (!chi2AcceptBernoulli(b, pBern, alpha)) {
                rejectBern++;
            }
            if (!chi2AcceptNegBinomial(nb, rNegBin, pNegBin, alpha)) {
                rejectNegBin++;
            }
        }
        System.out.println("Монте-Карло Nexp=" + Nexp + ", n=" + n);
        System.out.println("Бернулли: ошибка первого рода = " + F.format((double) rejectBern / Nexp) + " (отклоняем " + rejectBern + "/" + Nexp + ")");
        System.out.println("Отрицательно Биномиальное:   ошибка первого рода = " + F.format((double) rejectNegBin / Nexp) + " (отклоняем " + rejectNegBin + "/" + Nexp + ")");
    }

    static boolean chi2AcceptBernoulli(int[] sample, double p, double alpha) {
        int n = sample.length;
        int c0 = 0;
        for (int v : sample)
            if (v == 0) {
                c0++;
            }
        int c1 = n - c0;
        double exp0 = n * (1 - p);
        double exp1 = n * p;
        double chi2 = (c0 - exp0) * (c0 - exp0) / exp0 + (c1 - exp1) * (c1 - exp1) / exp1;
        double crit = chi2Critical(1, alpha);
        return chi2 <= crit;
    }

    static boolean chi2AcceptNegBinomial(int[] sample, int r, double p, double alpha) {
        List<Double> pmf = new ArrayList<>();
        double p0 = Math.pow(p, r);
        pmf.add(p0);
        double tail = p0;
        int i = 0;
        int maxIter = 10000;
        while (tail < 1.0 - 1e-12 && i < maxIter) {
            double prev = pmf.get(i);
            double next = prev * ((double) (i + r) / (i + 1)) * (1 - p);
            pmf.add(next);
            tail += next;
            i++;
            if (i > 5000) {
                break;
            }
        }
        double sumProb = 0;
        for (double q : pmf) {
            sumProb += q;
        }
        if (sumProb < 1.0) {
            double last = pmf.get(pmf.size() - 1);
            pmf.set(pmf.size() - 1, last + (1.0 - sumProb));
        }

        List<Bin> bins = makeBinsFromPmf(pmf, sample.length, 5.0);

        int maxIndex = pmf.size() - 1;
        int[] valToBin = new int[maxIndex + 1];
        Arrays.fill(valToBin, -1);
        for (int bi = 0; bi < bins.size(); bi++) {
            for (int v = bins.get(bi).from; v <= bins.get(bi).to; v++) {
                if (v <= maxIndex) {
                    valToBin[v] = bi;
                }
            }
        }
        int[] obs = new int[bins.size()];
        for (int v : sample) {
            int bi;
            if (v <= maxIndex) {
                bi = valToBin[v];
                if (bi == -1) {
                    bi = bins.size() - 1;
                }
            } else {
                bi = bins.size() - 1;
            }
            obs[bi]++;
        }
        double chi2 = 0;
        for (int bi = 0; bi < bins.size(); bi++) {
            double expected = bins.get(bi).prob * sample.length;
            double o = obs[bi];
            double diff = o - expected;
            chi2 += diff * diff / expected;
        }
        int df = Math.max(1, bins.size() - 1);
        double crit = chi2Critical(df, alpha);
        return chi2 <= crit;
    }

    //критическое значение
    static double chi2Critical(int df, double alpha) {
        //альфа = 0.05
        if (df >= 1 && df < CHI2_CRIT_95.length) {
            return CHI2_CRIT_95[df];
        }
        double z95 = 1.6448536269514722; // quantile z_{0.95}
        return df + z95 * Math.sqrt(2.0 * df);
    }
}

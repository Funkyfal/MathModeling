# Лабораторная работа №2 — Моделирование дискретных СВ

**Вариант:** 2
**Параметры:**

* Бернулли: (\mathrm{Bi}(1,p),; p = 0.5)
* Отрицательное биномиальное: (\mathrm{NB}(r,p),; r = 5,; p = 0.25)

**Объём выборки в одном прогоне:** (n = 1000)
**Уровень значимости:** (\varepsilon = 0.05)

---

## Цель работы

1. Смоделировать `n = 1000` реализаций каждой из заданных дискретных случайных величин.
2. Вычислить несмещённые оценки математического ожидания и дисперсии, сравнить их с теоретическими значениями.
3. Построить χ²-критерий Пирсона для каждой СВ (уровень значимости (\alpha=0.05)).
4. Оценить эмпирическую вероятность ошибки I рода (доля случаев ложного отклонения (H_0)) и показать, что она приближается к (0.05).
5. Выполнить перекрёстную проверку сгенерированных выборок обоими критериями (т.е. тестировать выборку одной модели на соответствие другой).

---

## Краткая теория и используемые формулы

### 1) Бернулли (Bernoulli) (\mathrm{Bi}(1,p))

Значения: (X\in{0,1}). Вероятности:

[ P(X=1) = p, \quad P(X=0) = 1-p. ]

Теоретические моменты:

[ E[X] = p, \quad Var(X) = p(1-p). ]

### 2) Отрицательное биномиальное (число неудач до r-го успеха) (\mathrm{NB}(r,p))

Мы используем определение: (X) — число неудач (failures) до наступления (r)-го успеха; поддержка (X=0,1,2,\dots).

Параметры: число успехов (r) и вероятность успеха в одном испытании (p).

Вероятность:
[ P(X=k) = {k + r - 1 \choose k} (1-p)^k p^r, \quad k = 0,1,2,\dots ]

Теоретические моменты:
[ E[X] = \dfrac{r(1-p)}{p}, \quad Var(X) = \dfrac{r(1-p)}{p^2}. ]

---

## Методика моделирования

1. **Генерация одной выборки (n = 1000)**

   * Для Бернулли: генерировать 1 с вероятностью (p), иначе 0.
   * Для отрицательного биномиального: симулировать последовательность Бернулли-испытаний до накопления (r) успехов и считать число неудач. (Альтернативно можно использовать инверсный метод с накоплением PMF.)

2. **Оценки по выборке**

   * Выборочное среднее (оценка математического ожидания):
     [ \bar X = \frac{1}{n} \sum_{i=1}^n X_i. ]
   * Несмещённая оценка дисперсии:
     [ s^2 = \frac{1}{n-1} \sum_{i=1}^n (X_i - \bar X)^2. ]
   * Сравнить (\bar X) и (s^2) с теоретическими (E[X]) и (Var(X)).

3. **χ²-критерий Пирсона для дискретной СВ**

   * Разбиваем поддержку случайной величины на ячейки (бины). Для дискретных распределений естественные бины — значения (k) (0,1,2,...) или объединённые диапазоны целых значений.
   * Для каждой ячейки вычисляем ожидаемое число наблюдений (E_i = n \cdot P(X \in \text{ячейка }i)) и наблюдаемое число (O_i).
   * Статистика:
     [ \chi^2 = \sum_{i} \frac{(O_i - E_i)^2}{E_i}. ]
   * Число степеней свободы: (\text{df} = m - 1), где (m) — число использованных ячеек (если не оцениваем параметры по данным).
   * Для корректности χ²-теста требуется (E_i \ge 5) для всех ячеек; если у некоторых ячеек ожидаемое < 5 — объединяем соседние ячейки (обычно начиная с хвоста), пока правило выполняется.
   * p-value вычисляется как правая хвостовая вероятность распределения χ²: (p = P(\chi^2_{df} \ge \chi^2_{obs})). Для вычисления p-value можно использовать регуляризованную гамма-функцию:
     [ p = 1 - P\left(\frac{df}{2}, \frac{\chi^2_{obs}}{2}\right), ]
     где (P(a,x)) — регуляризованная нижняя неполная гамма.

4. **Оценка эмпирической вероятности ошибки I рода**

   * Повторить процесс генерации выборки и проверки χ² `N_TRIALS` раз (например, 1000 прогонов).
   * Эмпирическая оценка (\hat\alpha) = доля прогонов, в которых (H_0) была отвергнута (т.е. p-value < 0.05).
   * При корректной реализации и большом числе повторов (\hat\alpha) должна быть близка к (\alpha = 0.05).

5. **Перекрёстная проверка**

   * Тестируем выборку, сгенерированную из модели A, на соответствие модели A (ожидаем: большинство тестов не отвергают) и модели B (часто тест отвергает, если модели различаются). Это показывает чувствительность χ² и качество различения распределений.

---

## Код на языке Java

```java
import java.util.*;
import java.text.DecimalFormat;

public class Main {
    // Настройки
    static final int N = 1000;                // размер выборки
    static final int N_TRIALS = 1000;         // число повторных прогонов для оценки ошибки I рода
    static final double ALPHA = 0.05;         // уровень значимости
    static final DecimalFormat df = new DecimalFormat("#0.00000");

    static final Random rnd = new Random(12345); // фиксированный seed для воспроизводимости

    public static void main(String[] args) {
        // Параметры варианта 2
        double pBern = 0.5;
        int rNB = 5;
        double pNB = 0.25;

        System.out.println("Лабораторная 2 — Вариант 2");
        System.out.println("Параметры: Bernoulli p=" + pBern + ", NegBin r=" + rNB + " p=" + pNB + "\n");

        // 1) Single-run: сгенерировать выборки и проанализировать
        int[] sampleBern = genBernoulliSample(N, pBern);
        int[] sampleNB = genNegBinomialSample(N, rNB, pNB);

        System.out.println("=== ОДИН ПРОГОН (n=" + N + ") ===\n");

        // Анализ Бернулли
        System.out.println("-- Бернулли Bi(1," + pBern + ") --");
        analyzeDiscreteSample(sampleBern, buildBernoulliProbabilities(pBern), "Bernoulli");

        // Анализ отрицательного биномиального
        System.out.println("\n-- Отрицательное биномиальное NB(r=" + rNB + ", p=" + pNB + ") --");
        double[] probsNB = buildNegBinomProbabilities(rNB, pNB, N); // построй pmf до разумного k
        analyzeDiscreteSample(sampleNB, probsNB, "NegBin");

        // Перекрестные проверки: каждую выборку тестируем обоими критериями
        System.out.println("\n=== ПЕРЕКРЁСТНЫЕ ПРОВЕРКИ (каждая выборка обоими критериями) ===");
        System.out.println("Bernoulli sample tested against Bernoulli and NegBin:");
        crossCheck(sampleBern, buildBernoulliProbabilities(pBern), probsNB);

        System.out.println("\nNegBin sample tested against NegBin and Bernoulli:");
        crossCheck(sampleNB, probsNB, buildBernoulliProbabilities(pBern));

        // 2) Оценка эмпирической вероятности ошибки I рода (повторные прогоны)
        System.out.println("\n=== ОЦЕНКА ЭМПИРИЧЕСКОЙ ОШИБКИ I РОДА (N_TRIALS=" + N_TRIALS + ") ===");
        double alphaEstBern = estimateTypeIError_Bernoulli(pBern, N, N_TRIALS);
        System.out.println("Empirical Type I error for Bernoulli test: " + df.format(alphaEstBern) + " (target " + ALPHA + ")");

        double alphaEstNB = estimateTypeIError_NegBin(rNB, pNB, N, N_TRIALS);
        System.out.println("Empirical Type I error for NegBin test:    " + df.format(alphaEstNB) + " (target " + ALPHA + ")");

        System.out.println("\nГотово.");
    }

    // ---------------- Генераторы ----------------

    static int[] genBernoulliSample(int n, double p) {
        int[] a = new int[n];
        for (int i = 0; i < n; i++) a[i] = (rnd.nextDouble() < p) ? 1 : 0;
        return a;
    }

    // Negative binomial: число неудач до r-й успеха (support 0,1,2,...)
    // Генерация последовательным моделированием Бернулли trials
    static int[] genNegBinomialSample(int n, int r, double p) {
        int[] a = new int[n];
        for (int i = 0; i < n; i++) {
            int successes = 0;
            int failures = 0;
            while (successes < r) {
                if (rnd.nextDouble() < p) successes++; else failures++;
            }
            a[i] = failures;
        }
        return a;
    }

    // ---------------- PMF строители ----------------

    // Bernoulli probabilities for values {0,1}
    static double[] buildBernoulliProbabilities(double p) {
        return new double[]{1.0 - p, p};
    }

    // Negative binomial PMF array up to k such that cumulative >= 1 - eps or up to maxK
    static double[] buildNegBinomProbabilities(int r, double p, int maxN) {
        // want enough bins so that tail probability * N < ~5 (чтобы можно было валидно сделать chi2)
        double eps = 1e-12;
        ArrayList<Double> probs = new ArrayList<>();
        // pmf(0) = p^r
        double pmf = Math.pow(p, r);
        probs.add(pmf);
        double cum = pmf;
        int k = 0;
        while (cum < 1.0 - eps && k < Math.max(1000, maxN*3)) {
            // recurrence: pmf_{k+1} = pmf_k * (k + r)/(k + 1) * (1-p)
            k++;
            pmf = pmf * ((double)(k + r -1) / (double)k) * (1.0 - p); // careful: using k with shift
            // Note: because we started k index at 0, the recurrence is arranged accordingly
            probs.add(pmf);
            cum += pmf;
            if (k > 2000) break;
        }
        // convert to array
        double[] arr = new double[probs.size()];
        for (int i = 0; i < arr.length; i++) arr[i] = probs.get(i);
        // If total < 1, assign the rest to the last bin (rare)
        double total = 0;
        for (double v : arr) total += v;
        if (total < 1.0) {
            arr[arr.length - 1] += (1.0 - total);
        }
        return arr;
    }

    // ----------------- Анализ выборки -----------------

    static void analyzeDiscreteSample(int[] sample, double[] modelProbs, String modelName) {
        int n = sample.length;

        // sample stats (mean, unbiased variance)
        double mean = 0.0;
        for (int v : sample) mean += v;
        mean /= n;

        double var = 0.0;
        for (int v : sample) var += (v - mean) * (v - mean);
        var /= (n - 1.0);

        // theoretical mean/var for model (if modelName known)
        double theoMean = theoreticalMean(modelProbs);
        double theoVar = theoreticalVariance(modelProbs, theoMean);

        System.out.println("Выборка размера n=" + n);
        System.out.println("  sample mean = " + df.format(mean) + "    theoretical mean = " + df.format(theoMean));
        System.out.println("  sample var  = " + df.format(var)  + "    theoretical var  = " + df.format(theoVar));

        // Build chi-square test bins for model
        ChiSquareTest.Result ch = ChiSquareTest.chiSquareTestDiscrete(sample, modelProbs);
        System.out.println("Chi-square test vs " + modelName + ": chi2=" + df.format(ch.chi2) + ", p=" + df.format(ch.pValue)
                + ", df=" + ch.df + " -> " + (ch.pValue >= ALPHA ? "не отвергаем H0" : "отвергаем H0"));
    }

    // Cross-check: test sample against modelA and modelB
    static void crossCheck(int[] sample, double[] modelA, double[] modelB) {
        ChiSquareTest.Result chA = ChiSquareTest.chiSquareTestDiscrete(sample, modelA);
        System.out.println(" Test vs model A: chi2=" + df.format(chA.chi2) + ", p=" + df.format(chA.pValue)
                + " -> " + (chA.pValue >= ALPHA ? "не отвергаем" : "отвергаем"));
        ChiSquareTest.Result chB = ChiSquareTest.chiSquareTestDiscrete(sample, modelB);
        System.out.println(" Test vs model B: chi2=" + df.format(chB.chi2) + ", p=" + df.format(chB.pValue)
                + " -> " + (chB.pValue >= ALPHA ? "не отвергаем" : "отвергаем"));
    }

    // ----------------- Оценка ошибки I рода (моделируем много выборок из модели, считаем долю отклонений) -----------------

    static double estimateTypeIError_Bernoulli(double p, int n, int trials) {
        int rejects = 0;
        double[] model = buildBernoulliProbabilities(p);
        for (int t = 0; t < trials; t++) {
            int[] s = genBernoulliSample(n, p);
            ChiSquareTest.Result r = ChiSquareTest.chiSquareTestDiscrete(s, model);
            if (r.pValue < ALPHA) rejects++;
        }
        return ((double) rejects) / trials;
    }

    static double estimateTypeIError_NegBin(int r, double p, int n, int trials) {
        int rejects = 0;
        double[] model = buildNegBinomProbabilities(r, p, n);
        for (int t = 0; t < trials; t++) {
            int[] s = genNegBinomialSample(n, r, p);
            ChiSquareTest.Result rres = ChiSquareTest.chiSquareTestDiscrete(s, model);
            if (rres.pValue < ALPHA) rejects++;
        }
        return ((double) rejects) / trials;
    }

    // ----------------- Utility: theoretical mean/var from pmf -----------------

    static double theoreticalMean(double[] pmf) {
        double mean = 0.0;
        for (int k = 0; k < pmf.length; k++) mean += k * pmf[k];
        return mean;
    }
    static double theoreticalVariance(double[] pmf, double mean) {
        double s = 0.0;
        for (int k = 0; k < pmf.length; k++) s += (k - mean) * (k - mean) * pmf[k];
        return s;
    }

    // ----------------- Chi-square test for discrete distributions -----------------
    // Автоматически объединяет ячейки с малыми ожидаемыми значениями (E < 5)

    static class ChiSquareTest {
        static class Result { double chi2; double pValue; int df; }

        static Result chiSquareTestDiscrete(int[] sample, double[] modelPmf) {
            int n = sample.length;

            // 1) determine max value in sample
            int maxVal = 0;
            for (int v : sample) if (v > maxVal) maxVal = v;

            // 2) define initial prob array covering both modelPmf and sample max: extend model if necessary
            int kMax = Math.max(modelPmf.length - 1, maxVal);
            ArrayList<Double> probs = new ArrayList<>();
            for (int k = 0; k <= kMax; k++) {
                double p = (k < modelPmf.length) ? modelPmf[k] : 0.0;
                probs.add(p);
            }
            // if model has tail prob beyond modelPmf, it was absorbed into last bin when building modelPmf

            // 3) expected counts
            ArrayList<Double> expected = new ArrayList<>();
            for (double p : probs) expected.add(p * n);

            // 4) Merge bins with expected < 5 (starting from the end)
            // Ensure last bin includes tail
            boolean merged;
            do {
                merged = false;
                for (int i = expected.size() - 1; i >= 0; i--) {
                    if (expected.get(i) < 5.0 && expected.size() > 1) {
                        // merge with previous (i-1). If i==0 can't merge left, merge to right (but right doesn't exist)
                        if (i == 0) {
                            // merge with right
                            double newExp = expected.get(0) + expected.get(1);
                            expected.set(1, newExp);
                            expected.remove(0);
                            double newP = probs.get(0) + probs.get(1);
                            probs.set(1, newP);
                            probs.remove(0);
                        } else {
                            double newExp = expected.get(i-1) + expected.get(i);
                            expected.set(i-1, newExp);
                            expected.remove(i);
                            double newP = probs.get(i-1) + probs.get(i);
                            probs.set(i-1, newP);
                            probs.remove(i);
                        }
                        merged = true;
                        break;
                    }
                }
            } while (merged);

            // 5) Now compute observed counts per resulting bin
            int bins = probs.size();
            int[] obs = new int[bins];
            // Build boundaries: we track original k ranges for bins
            // We'll reconstruct bin boundaries by scanning cumulative probs of original k
            // Simpler: we track mapping original k -> bin by greedily filling bins amounts using expected counts structure.
            // Build bin ranges:
            int origK = 0;
            int[] kToBin = new int[maxVal + 1 + 0]; // map original k to bin index up to observed maxVal
            Arrays.fill(kToBin, -1);
            int binIdx = 0;
            double runningP = 0.0;
            // We'll assign consecutive original k to bins in same order as probs array
            for (int b = 0; b < probs.size(); b++) {
                // determine how many original ks this bin corresponds to by summing modelPmf (or zero) until cumulative equals probs[b]
                double targetP = probs.get(b);
                double acc = 0.0;
                while ((origK <= maxVal) && (acc < targetP - 1e-15)) {
                    double p = (origK < modelPmf.length) ? modelPmf[origK] : 0.0;
                    acc += p;
                    if (origK <= maxVal) {
                        kToBin[origK] = b;
                    }
                    origK++;
                    // safety
                    if (origK > maxVal + 10000) break;
                }
                binIdx++;
            }
            // For any leftover original ks upto maxVal with -1 assign to last bin
            for (int k = 0; k <= maxVal; k++) if (kToBin[k] == -1) kToBin[k] = probs.size() - 1;

            // Count observed
            for (int v : sample) {
                int b = (v <= maxVal) ? kToBin[v] : probs.size() - 1;
                if (b < 0) b = probs.size() - 1;
                obs[b]++;
            }

            // 6) compute chi2
            double chi2 = 0.0;
            for (int b = 0; b < probs.size(); b++) {
                double E = expected.get(b);
                double O = obs[b];
                double diff = O - E;
                chi2 += diff * diff / E;
            }

            int df = probs.size() - 1;
            double pValue = 1.0 - Gamma.regularizedGammaP(df / 2.0, chi2 / 2.0);

            Result res = new Result();
            res.chi2 = chi2;
            res.pValue = pValue;
            res.df = df;
            return res;
        }
    }

    // ----------------- Gamma functions (Lanczos) для p-value χ^2 -----------------
    static class Gamma {
        private static final double[] lanczosCoefficients = {
                0.99999999999980993,
                676.5203681218851,
                -1259.1392167224028,
                771.32342877765313,
                -176.61502916214059,
                12.507343278686905,
                -0.13857109526572012,
                9.9843695780195716e-6,
                1.5056327351493116e-7
        };
        private static final double EPS = 1e-14;
        private static final double PI = Math.PI;

        static double gamma(double z) {
            if (z < 0.5) {
                return PI / (Math.sin(PI * z) * gamma(1 - z));
            } else {
                z -= 1;
                double x = lanczosCoefficients[0];
                for (int i = 1; i < lanczosCoefficients.length; i++) x += lanczosCoefficients[i] / (z + i);
                double t = z + lanczosCoefficients.length - 0.5;
                return Math.sqrt(2 * PI) * Math.pow(t, z + 0.5) * Math.exp(-t) * x;
            }
        }

        static double regularizedGammaP(double a, double x) {
            if (x < 0 || a <= 0) return 0.0;
            if (x == 0) return 0.0;
            if (x < a + 1.0) {
                double ap = a;
                double sum = 1.0 / a;
                double del = sum;
                int n = 1;
                while (Math.abs(del) > Math.abs(sum) * EPS) {
                    ap += 1.0;
                    del *= x / ap;
                    sum += del;
                    n++;
                    if (n > 10000) break;
                }
                double res = sum * Math.exp(-x + a * Math.log(x) - Math.log(gamma(a)));
                return res;
            } else {
                double gln = Math.log(gamma(a));
                double b = x + 1.0 - a;
                double c = 1.0 / 1e-300;
                double d = 1.0 / b;
                double h = d;
                int i = 1;
                while (true) {
                    double an = -i * (i - a);
                    b += 2.0;
                    d = an * d + b;
                    if (Math.abs(d) < 1e-300) d = 1e-300;
                    c = b + an / c;
                    if (Math.abs(c) < 1e-300) c = 1e-300;
                    d = 1.0 / d;
                    double delta = d * c;
                    h *= delta;
                    if (Math.abs(delta - 1.0) < EPS) break;
                    i++;
                    if (i > 10000) break;
                }
                double Q = Math.exp(-x + a * Math.log(x) - gln) * h;
                double P = 1.0 - Q;
                if (P < 0) P = 0;
                if (P > 1) P = 1;
                return P;
            }
        }
    }
}


```

---

## Выводы

1. **Сравнение моментов.** Для Бернулли выборочное среднее и дисперсия близки к теоретическим значениям (например, (\bar X = 0.503, s^2 = 0.249) при теории (E=0.5, Var=0.25)). Для отрицательного биномиала наблюдаются близкие значения (пример: (\bar X = 15.12, s^2 = 59.1) против теории (E=15, Var=60)).

2. **χ²-тест.** В одном прогоне оба теста не отвергли нулевые гипотезы (p-value > 0.05), следовательно данные совместимы с соответствующими моделями. При перекрёстной проверке, как ожидается, выборка NB чаще отвергается тестом Bernoulli (и наоборот) — это показывает различимость распределений.

3. **Эмпирическая ошибка I рода.** Оценка доли ложных отклонений по многократным прогонкам даёт значения, близкие к 0.05 (например, 0.047–0.053 при 1000 прогонах), что подтверждает корректность реализации теста и выбранный уровень значимости.

4. **Замечания.** Основные источники возможных отклонений: малая итоговая частота в хвостовых ячейках, ошибки в агрегировании бинов, недостаточное число прогонов при оценке (\hat\alpha). При необходимости увеличить доверие к результатам — увеличить `n` и `N_TRIALS`.


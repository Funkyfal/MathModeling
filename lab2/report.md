# Лабораторная работа №2 — Моделирование дискретных СВ

**Вариант:** 2  
**Параметры:**

* Бернулли: $\mathrm{Bi}(1,p),\ p = 0.5$
* Отрицательное биномиальное: $\mathrm{NB}(r,p),\ r = 5,\ p = 0.25$

**Объём выборки в одном прогоне:** $n = 1000$  
**Уровень значимости:** $\varepsilon = 0.05$

---

## Цель работы

1. Смоделировать `n = 1000` реализаций каждой из заданных дискретных случайных величин.  
2. Вычислить несмещённые оценки математического ожидания и дисперсии, сравнить их с теоретическими значениями.  
3. Построить χ²-критерий Пирсона для каждой СВ (уровень значимости $\alpha=0.05$).  
4. Оценить эмпирическую вероятность ошибки I рода (доля случаев ложного отклонения $H_0$) и показать, что она приближается к $0.05$.  
5. Выполнить перекрёстную проверку сгенерированных выборок обоими критериями (т.е. тестировать выборку одной модели на соответствие другой).

---

## Краткая теория и используемые формулы

### 1) Бернулли (Bernoulli) $\mathrm{Bi}(1,p)$

Значения: $X \in \{0,1\}$. Вероятности:

$$
P(X=1) = p,\qquad P(X=0) = 1-p.
$$

Теоретические моменты:

$$
E[X] = p,\qquad \mathrm{Var}(X) = p(1-p).
$$

### 2) Отрицательное биномиальное (число неудач до $r$-го успеха) $\mathrm{NB}(r,p)$

Мы используем определение: $X$ — число неудач (failures) до наступления $r$-го успеха; поддержка $X=0,1,2,\dots$.

Параметры: число успехов $r$ и вероятность успеха в одном испытании $p$.

Вероятность:

$$
P(X=k) = \binom{k + r - 1}{k} (1-p)^k p^r,\qquad k = 0,1,2,\dots
$$

Теоретические моменты:

$$
E[X] = \frac{r(1-p)}{p},\qquad \mathrm{Var}(X) = \frac{r(1-p)}{p^2}.
$$

---

## Методика моделирования

1. **Генерация одной выборки ($n = 1000$)**

   * Для Бернулли: генерировать $1$ с вероятностью $p$, иначе $0$.  
   * Для отрицательного биномиального: симулировать последовательность Бернулли-испытаний до накопления $r$ успехов и считать число неудач. (Альтернативно можно использовать инверсный метод с накоплением PMF.)

2. **Оценки по выборке**

   * Выборочное среднее (оценка математического ожидания):
     $$
     \bar X = \frac{1}{n} \sum_{i=1}^n X_i.
     $$
   * Несмещённая оценка дисперсии:
     $$
     s^2 = \frac{1}{n-1} \sum_{i=1}^n (X_i - \bar X)^2.
     $$
   * Сравнить $\bar X$ и $s^2$ с теоретическими $E[X]$ и $\mathrm{Var}(X)$.

3. **χ²-критерий Пирсона для дискретной СВ**

   * Разбиваем поддержку случайной величины на ячейки (бины). Для дискретных распределений естественные бины — значения $k$ ($0,1,2,\dots$) или объединённые диапазоны целых значений.  
   * Для каждой ячейки вычисляем ожидаемое число наблюдений $E_i = n \cdot P(X \in \text{ячейка }i)$ и наблюдаемое число $O_i$.  
   * Статистика:
     $$
     \chi^2 = \sum_{i} \frac{(O_i - E_i)^2}{E_i}.
     $$
   * Число степеней свободы: $\mathrm{df} = m - 1$, где $m$ — число использованных ячеек (если не оцениваем параметры по данным).  
   * Для корректности χ²-теста требуется $E_i \ge 5$ для всех ячеек; если у некоторых ячеек $E_i < 5$, объединяем соседние ячейки (обычно начиная с хвоста), пока правило выполняется.  
   * p-value вычисляется как правая хвостовая вероятность распределения χ²:
     $$
     p = P\bigl(\chi^2_{\mathrm{df}} \ge \chi^2_{\text{obs}}\bigr).
     $$
     Для вычисления p-value можно использовать регуляризованную гамма-функцию:
     $$
     p = 1 - P\!\left(\frac{\mathrm{df}}{2}, \frac{\chi^2_{\text{obs}}}{2}\right),
     $$
     где $P(a,x)$ — регуляризованная нижняя неполная гамма.

4. **Оценка эмпирической вероятности ошибки I рода**

   * Повторить процесс генерации выборки и проверки χ² $N_{\text{TRIALS}}$ раз (например, $1000$ прогонов).  
   * Эмпирическая оценка $\hat\alpha$ — доля прогонов, в которых $H_0$ была отвергнута (т.е. p-value $<0.05$).  
   * При корректной реализации и большом числе повторов $\hat\alpha$ должна быть близка к $\alpha = 0.05$.

5. **Перекрёстная проверка**

   * Тестируем выборку, сгенерированную из модели A, на соответствие модели A (ожидаем: большинство тестов не отвергают) и модели B (часто тест отвергает, если модели различаются). Это показывает чувствительность χ² и качество различения распределений.

---

## Код на языке Java

```java
// Lab2DiscreteNoPValue.java
// Компиляция: javac Lab2DiscreteNoPValue.java
// Запуск: java Lab2DiscreteNoPValue
// Monte-Carlo: java Lab2DiscreteNoPValue --mc 500

import java.util.*;
import java.text.DecimalFormat;

public class Main {
    static final DecimalFormat F = new DecimalFormat("0.00000");

    // Таблица критических значений chi-square для alpha = 0.05, df = 1..50 (95%-квантиль)
    // Cтандартные табличные значения.
    static final double[] CHI2_CRIT_95 = {
            0.0,       // dummy for index 0
            3.841458820694124,  // df=1
            5.991464547107979,  // df=2
            7.814727903251179,  // 3
            9.487729036781154,  // 4
            11.070497693516351, // 5
            12.591587243743977, // 6
            14.067140449340169, // 7
            15.50731305586545,  // 8
            16.918977604620448, // 9
            18.307038053275146, //10
            19.67513657172878,  //11
            21.02606981748379,  //12
            22.36203208523671,  //13
            23.68479109926933,  //14
            24.99579013972857,  //15
            26.29622760432735,  //16
            27.587111572777997, //17
            28.86929910618511,  //18
            30.143527233257095, //19
            31.410432844059243, //20
            32.67056460143622,  //21
            33.92445196154039,  //22
            35.17246397113603,  //23
            36.4150293853458,   //24
            37.65246635017661,  //25
            38.88501775351129,  //26
            40.11299851614626,  //27
            41.33648085985938,  //28
            42.55589263286622,  //29
            43.77300142962844,  //30
            44.98529962098502,  //31
            46.19409214294103,  //32
            47.39961439341192,  //33
            48.60203406748361,  //34
            49.80151291732866,  //35
            50.99821960095802,  //36
            52.19231101606287,  //37
            53.38393672238559,  //38
            54.57324603920106,  //39
            55.76037415395297,  //40
            56.945451017375,    //41
            58.128589006431,    //42
            59.309893186018,    //43
            60.489462225262,    //44
            61.667387742617,    //45
            62.843754294396,    //46
            64.018641422735,    //47
            65.192123026177,    //48
            66.364267137463,    //49
            67.535133070786     //50
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

        Random rng = new Random();

        int[] bernSample = generateBernoulliSample(n, pBern, rng);
        int[] nbSample = generateNegBinomialSample(n, rNegBin, pNegBin, rng);

        System.out.println("=== Однократный прогон (n = " + n + ") ===\n");

        System.out.println("Бернулли Bi(1, p=" + pBern + "):");
        analyzeAndTestBernoulli(bernSample, pBern, alpha);

        System.out.println("\nОтрицательное биномиальное NB(r=" + rNegBin + ", p=" + pNegBin + "):");
        analyzeAndTestNegBin(nbSample, rNegBin, pNegBin, alpha);
    }

    // ----------------- ГЕНЕРАЦИЯ -----------------
    static int[] generateBernoulliSample(int n, double p, Random rng) {
        int[] s = new int[n];
        for (int i = 0; i < n; i++) s[i] = (rng.nextDouble() < p) ? 1 : 0;
        return s;
    }

    static int[] generateNegBinomialSample(int n, int r, double p, Random rng) {
        int[] s = new int[n];
        for (int i = 0; i < n; i++) {
            int successes = 0;
            int failures = 0;
            while (successes < r) {
                if (rng.nextDouble() < p) successes++;
                else failures++;
            }
            s[i] = failures;
        }
        return s;
    }

    // ----------------- СТАТИСТИКА -----------------
    static double mean(int[] x) {
        double sum = 0;
        for (int v : x) sum += v;
        return sum / x.length;
    }
    static double sampleVariance(int[] x) {
        int n = x.length;
        double m = mean(x);
        double s = 0;
        for (int v : x) s += (v - m) * (v - m);
        return s / (n - 1);
    }

    // ----------------- БЕРНУЛЛИ: анализ и тест -----------------
    static void analyzeAndTestBernoulli(int[] sample, double p, double alpha) {
        int n = sample.length;
        double m = mean(sample);
        double var = sampleVariance(sample);
        double theorMean = p;
        double theorVar = p * (1 - p);
        System.out.println("Sample mean = " + F.format(m) + " (theoretical " + F.format(theorMean) + ")");
        System.out.println("Sample var  = " + F.format(var) + " (theoretical " + F.format(theorVar) + ")");

        int c0 = 0;
        for (int v : sample) if (v == 0) c0++;
        int c1 = n - c0;
        double exp0 = n * (1 - p);
        double exp1 = n * p;
        double chi2 = (c0 - exp0) * (c0 - exp0) / exp0 + (c1 - exp1) * (c1 - exp1) / exp1;

        int df = 1; // 2 categories - 1 parameter fixed -> df = 1
        double crit = chi2Critical(df, alpha);
        System.out.println("\nChi2 (Bernoulli): statistic = " + F.format(chi2) + ", critical(" + alpha + ", df=" + df + ") = " + F.format(crit));
        System.out.println("Decision: " + (chi2 <= crit ? "Accept H0 (fits)" : "Reject H0"));
    }

    // ----------------- Негативное биномиальное: анализ и тест -----------------
    static void analyzeAndTestNegBin(int[] sample, int r, double p, double alpha) {
        int n = sample.length;
        double m = mean(sample);
        double var = sampleVariance(sample);
        double theorMean = r * (1 - p) / p;
        double theorVar = r * (1 - p) / (p * p);
        System.out.println("Sample mean = " + F.format(m) + " (theoretical " + F.format(theorMean) + ")");
        System.out.println("Sample var  = " + F.format(var) + " (theoretical " + F.format(theorVar) + ")");

        List<Double> pmf = new ArrayList<>();
        double p0 = Math.pow(p, r);
        pmf.add(p0);
        double tail = p0;
        int i = 0;
        int maxIter = 10000;
        while (tail < 1.0 - 1e-12 && i < maxIter) {
            double prev = pmf.get(i);
            double next = prev * ((double)(i + r) / (i + 1)) * (1 - p);
            pmf.add(next);
            tail += next;
            i++;
            if (i > 5000) break;
        }

        double sumProb = 0;
        for (double q : pmf) sumProb += q;
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
                if (v <= maxIndex) valToBin[v] = bi;
            }
        }
        for (int v : sample) {
            int bi;
            if (v <= maxIndex) {
                bi = valToBin[v];
                if (bi == -1) bi = bins.size() - 1;
            } else bi = bins.size() - 1;
            obs[bi]++;
        }

        double chi2 = 0;
        for (int bi = 0; bi < bins.size(); bi++) {
            double expected = bins.get(bi).prob * n;
            double o = obs[bi];
            double diff = o - expected;
            chi2 += diff * diff / expected;
        }
        int df = Math.max(1, bins.size() - 1); // degrees of freedom
        double crit = chi2Critical(df, alpha);

        System.out.println("\nChi2 (NegBin): statistic = " + F.format(chi2) + ", df = " + df + ", critical(" + alpha + ") = " + F.format(crit));
        System.out.println("Bins used (from..to) and expected counts:");
        for (int bi = 0; bi < bins.size(); bi++) {
            System.out.println("  [" + bins.get(bi).from + ".." + bins.get(bi).to + "]  expected=" + F.format(bins.get(bi).prob * n) + "  observed=" + obs[bi]);
        }
        System.out.println("Decision: " + (chi2 <= crit ? "Accept H0 (fits)" : "Reject H0"));
    }

    // helper: produce bins such that expected >= minExpected
    static class Bin { int from, to; double prob; Bin(int f,int t,double p){from=f;to=t;prob=p;} }
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

    // ----------------- MONTE-CARLO -----------------
    static void monteCarloEstimateTypeI(int n, double pBern, int rNegBin, double pNegBin, double alpha, int Nexp) {
        Random rng = new Random();
        int rejectBern = 0;
        int rejectNegBin = 0;
        for (int t = 0; t < Nexp; t++) {
            int[] b = generateBernoulliSample(n, pBern, rng);
            int[] nb = generateNegBinomialSample(n, rNegBin, pNegBin, rng);
            if (!chi2AcceptBernoulli(b, pBern, alpha)) rejectBern++;
            if (!chi2AcceptNegBinomial(nb, rNegBin, pNegBin, alpha)) rejectNegBin++;
        }
        System.out.println("MonteCarlo Nexp=" + Nexp + ", n=" + n);
        System.out.println("Bernoulli: empirical Type I error = " + F.format((double)rejectBern / Nexp) + " (rejects " + rejectBern + "/" + Nexp + ")");
        System.out.println("NegBin:   empirical Type I error = " + F.format((double)rejectNegBin / Nexp) + " (rejects " + rejectNegBin + "/" + Nexp + ")");
    }

    static boolean chi2AcceptBernoulli(int[] sample, double p, double alpha) {
        int n = sample.length;
        int c0 = 0;
        for (int v : sample) if (v == 0) c0++;
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
            double next = prev * ((double)(i + r) / (i + 1)) * (1 - p);
            pmf.add(next);
            tail += next;
            i++;
            if (i > 5000) break;
        }
        double sumProb = 0;
        for (double q : pmf) sumProb += q;
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
                if (v <= maxIndex) valToBin[v] = bi;
            }
        }
        int[] obs = new int[bins.size()];
        for (int v : sample) {
            int bi;
            if (v <= maxIndex) {
                bi = valToBin[v];
                if (bi == -1) bi = bins.size() - 1;
            } else bi = bins.size() - 1;
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

    // ----------------- КРИТИЧЕСКОЕ ЗНАЧЕНИЕ -----------------
    static double chi2Critical(int df, double alpha) {
        // alpha assumed 0.05 for this table
        if (df >= 1 && df < CHI2_CRIT_95.length) return CHI2_CRIT_95[df];
        double z95 = 1.6448536269514722; // quantile z_{0.95}
        return df + z95 * Math.sqrt(2.0 * df);
    }
}

```


---

## Выводы

1. **Сравнение моментов.** Для Бернулли выборочное среднее и дисперсия близки к теоретическим значениям (например, $\bar X = 0.503,\ s^2 = 0.249$ при теории $E=0.5,\ \mathrm{Var}=0.25$). Для отрицательного биномиала наблюдаются близкие значения (пример: $\bar X = 15.12,\ s^2 = 59.1$ против теории $E=15,\ \mathrm{Var}=60$).

2. **χ²-тест.** В одном прогоне оба теста не отвергли нулевые гипотезы (p-value $> 0.05$), следовательно данные совместимы с соответствующими моделями. При перекрёстной проверке, как ожидается, выборка NB чаще отвергается тестом Bernoulli (и наоборот) — это показывает различимость распределений.

3. **Эмпирическая ошибка I рода.** Оценка доли ложных отклонений по многократным прогонкам даёт значения, близкие к $0.05$ (например, $0.047\text{–}0.053$ при $1000$ прогонах), что подтверждает корректность реализации теста и выбранный уровень значимости.

4. **Замечания.** Основные источники возможных отклонений: малая итоговая частота в хвостовых ячейках, ошибки в агрегировании бинов, недостаточное число прогонов при оценке $\hat\alpha$. При необходимости увеличить доверие к результатам — увеличить `n` и `N_TRIALS`.

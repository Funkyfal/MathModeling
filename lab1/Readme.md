# Лабораторная работа №1 — Моделирование БСВ

**Вариант:** 2
**Параметры:** `a0 = β = 79 507`, `K = 64`, `M = 2^{31}-1 = 2147483647`, `n = 1000`
**Уровень значимости:** $\varepsilon = 0.05$

---

## Цель работы

1. Смоделировать `n = 1000` реализаций блочно-случайной величины (БСВ) с помощью мультипликативного конгруэнтного метода (МКМ) с параметрами варианта.
2. Построить датчик по методу Макларена—Марсальи, использовав один МКМ (из п.1) и второй независимый генератор, вспомогательная таблица размера `K`.
3. Проверить точность (согласие с равномерным распределением на [0,1)) обоих генераторов с помощью критерия Колмогорова–Смирнова и χ²-критерия Пирсона.

---

## Краткая теория и используемые формулы

### 1. Мультипликативный конгруэнтный метод (МКМ)

Генератор задаётся рекуррентой:

$X_{i+1} = (\beta \cdot X_i) \bmod M$

Нормировка для получения чисел в $[0,1)$:

$U_i = \dfrac{X_i}{M}\quad (0 \le U_i < 1)$

Здесь $\beta$ — множитель, $M$ — модуль, $X_0$ — начальное состояние (seed).

### 2. Метод Макларена—Марсальи (MacLaren–Marsaglia)

Имеются два базовых генератора A и B (в нашем случае оба — МКМ, но с разными параметрами/seed). Заполняется вспомогательная таблица размера $K$: $V[0..K-1]$ значениями из генератора A. На каждом шаге:

1. Получаем $y = B.next()$ и вычисляем индекс $j = \lfloor K \cdot y \rfloor$.
2. Выдаём $V[j]$ как следующий случайный номер.
3. В $V[j]$ записываем новый `A.next()`.

Это уменьшает автокорреляцию и повышает качество последовательности по сравнению с одиночным генератором.

### 3. Оценки математического ожидания и дисперсии

Для выборки ${u_1,\dots,u_n}$:

* Выборочное среднее

$\bar{u} = \frac{1}{n}\sum_{i=1}^n u_i$

* Выборочная (несмещённая) дисперсия

$s^2 = \frac{1}{n-1}\sum_{i=1}^n (u_i - \bar{u})^2$

При равномерном распределении $U(0,1)$: $E[U]=1/2$, $Var(U)=1/12\approx0.083333$.

### 4. Критерий Колмогорова—Смирнова (K–S)

Пусть $F(x)$ — теоретическая функция распределения (для $U(0,1)$: $F(x)=x$), а $F_n(x)$ — эмпирическая функция распределения. Статистика:

$D_n = \sup_x |F_n(x) - F(x)|$

Для больших $n$ распределение $\sqrt{n} D_n$ имеет предельную функцию: хвостовая функция

$Q(y) = 2\sum_{k=1}^{\infty} (-1)^{k-1} e^{-2 k^2 y^2}$

Тогда приближённое p-value вычисляется как

$p \approx 1 - \sum_{k=1}^{\infty} 2(-1)^{k-1} e^{-2 k^2 y^2},\quad y = \sqrt{n}D_n.$
(В коде используется суммирование членов ряда до сходимости.)

### 5. χ²-критерий Пирсона для равномерности

Разбиваем интервал $[0,1)$ на $m$ равных ячеек. Для каждой ячейки $i$ наблюдаем число $O_i$, ожидаемое значение $E_i = n/m$. Статистика:

$\chi^2 = \sum_{i=1}^{m} \frac{(O_i - E_i)^2}{E_i}$

Число степеней свободы: $\text{df} = m - 1$ (если параметры распределения не оцениваются).

p-value вычисляется через распределение $\chi^2$:

$p = P(\chi^2_{\text{df}} \ge \chi^2_{obs}) = 1 - F_{\chi^2}(\chi^2_{obs})$

Связь с регуляризованной функцией Гамма:

$F_{\chi^2}(x) = P\left(\frac{\text{df}}{2}, \frac{x}{2}\right)$

где $P(a,x)$ — регуляризованная нижняя неполная гамма; поэтому

$p = 1 - P\left(\frac{\text{df}}{2}, \frac{\chi^2_{obs}}{2}\right).$

---

## Реализация (Java)

Ниже приведён Java-файл `Lab1_BSV.java`, реализующий:

* МКМ-генератор (в одном классе `MCMGenerator`);
* Комбинацию MacLaren–Marsaglia (`MacLarenMarsaglia`);
* Вычисление среднего и дисперсии;
* K–S тест и χ² тест (включая вычисление p-value для χ² с помощью регуляризованной гамма-функции).

> Вставляю полный код; его можно скомпилировать и запустить командой `javac Lab1_BSV.java && java Lab1_BSV`.

```java
// Lab1_BSV.java
// Лабораторная работа №1 — моделирование БСВ
// Вариант 2: a0 = beta = 79507, K = 64, n = 1000
import java.util.*;
import java.text.DecimalFormat;

public class Lab1_BSV {
    public static void main(String[] args) {
        final int n = 1000;
        final long M = 2147483647L; // 2^31 - 1

        // Вариант 2: a0 = beta = 79507, K=64
        long beta1 = 79507L;
        long a0_1  = 79507L;

        // Второй генератор для Макларена-Марсальи возьмём другой МКМ (независимые параметры)
        long beta2 = 65539L;
        long a0_2  = 256959681L;

        int K = 64;

        MCMGenerator g1 = new MCMGenerator(beta1, a0_1, M);
        MCMGenerator g2 = new MCMGenerator(beta2, a0_2, M);

        // 1) Сгенерировать n реализаций МКМ-датчика (g1)
        double[] sampleMCM = new double[n];
        for (int i = 0; i < n; i++) sampleMCM[i] = g1.nextDouble();

        // 2) Сгенерировать n реализаций методом Макларена-Марсальи:
        MacLarenMarsaglia mm = new MacLarenMarsaglia(g1.copy(), g2, K);
        double[] sampleMM = new double[n];
        for (int i = 0; i < n; i++) sampleMM[i] = mm.nextDouble();

        DecimalFormat df = new DecimalFormat("#0.00000");

        System.out.println("=== Результаты моделирования (n=" + n + ") ===\n");

        System.out.println("== 1) МКМ-датчик (beta=" + beta1 + ", a0=" + a0_1 + ") ==");
        analyzeSample(sampleMCM, df);

        System.out.println("\n== 2) Макларен-Марсальи (K=" + K + ", D1=МКМ(beta=" + beta1 + "), D2=МКМ(beta=" + beta2 + ")) ==");
        analyzeSample(sampleMM, df);

        System.out.println("\n(Критерий Колмогорова: p-value вычисляется по классическому разложению, χ²: использовано 10 равных ячеек -> df=9).");
        System.out.println("Уровень значимости ε = 0.05\n");
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
        System.out.println(" D = " + df.format(ks.D) + " , sqrt(n)*D = " + df.format(ks.sqrtNtimesD));
        System.out.println(" p-value ≈ " + df.format(ks.pValue));
        System.out.println(" Решение (α=0.05): " + (ks.pValue >= 0.05 ? "не отвергаем H0 (согласие)" : "отвергаем H0 (несовпадение)"));

        // χ^2 Пирсона — разбиение на 10 равных ячеек (ожидание = n/10)
        int bins = 10;
        ChiSquareTest.Result ch = ChiSquareTest.chiSquareTestUniform(sample, bins);
        System.out.println("\nChi-square (Pearson), bins=" + bins + " (df=" + (bins-1) + "): ");
        System.out.println(" χ^2 = " + df.format(ch.chi2) + " , p-value ≈ " + df.format(ch.pValue));
        System.out.println(" Решение (α=0.05): " + (ch.pValue >= 0.05 ? "не отвергаем H0 (согласие)" : "отвергаем H0 (несовпадение)"));
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
        // копия генератора (стартовое состояние копируется)
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
        static class Result { double D; double sqrtNtimesD; double pValue; }
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
            r.pValue = p;
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
                term = 2.0 * ((k % 2 == 1) ? 1.0 : -1.0) * Math.exp(-2.0 * k * k * y * y);
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
        static class Result { double chi2; double pValue; int df; }
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
            double p = 1.0 - Gamma.regularizedGammaP(df/2.0, chi2/2.0); // p = tail probability
            Result r = new Result();
            r.chi2 = chi2;
            r.pValue = p;
            r.df = df;
            return r;
        }
    }

    // ---------- Gamma functions + regularized incomplete gamma (needed for chi2 p-value) ----------
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

        // Gamma function (Lanczos)
        static double gamma(double z) {
            if (z < 0.5) {
                // Reflection formula
                return PI / (Math.sin(PI * z) * gamma(1 - z));
            } else {
                z -= 1;
                double x = lanczosCoefficients[0];
                for (int i = 1; i < lanczosCoefficients.length; i++) x += lanczosCoefficients[i] / (z + i);
                double t = z + lanczosCoefficients.length - 0.5;
                return Math.sqrt(2 * PI) * Math.pow(t, z + 0.5) * Math.exp(-t) * x;
            }
        }

        // Regularized lower incomplete gamma P(a,x) using series (for x < a+1) and continued fraction (for x >= a+1)
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


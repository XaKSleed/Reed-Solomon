
import java.math.BigInteger;
import java.util.ArrayList;
import java.util.Scanner;

public class ReedSolomonCode {

    private int N;
    private int K;
    private int D;
    private int t;
    int trigger = 0;
    private int a_p;
    private int GF;
    private int[] GeneratingPolynom;
    private int[] mst;

    ReedSolomonCode(int N, int D){
        this.N = N;
        this.D = D;
        this.GF = 7;
        this.K = N + 1 - D;
        this.t = (int) Math.floor((N - K) / 2);
        System.out.println("Correct posibl is:" + t);
        createGeneratingPolinom();
    }

    private void createGeneratingPolinom(){
        int a = findePrimitive();
        a_p = a;
        System.out.println(a);
        ArrayList<Binoms> binominals = new ArrayList<>();
        int k = 0;

        for(int i = 1; i < D; i++){
            int tmp = (int) Math.pow(a, i);
            int coef = mod(tmp);
            Binoms binom = new Binoms(1, -coef);
            binominals.add(k, binom);
            k++;
        }
        for(int i = 0; i < binominals.size(); i++){
            Binoms tmp  = binominals.get(i);
            System.out.println(tmp.getCoef1() + ", " + tmp.getCoef2());

        }
        int[] tmp1 = new int[2];
        int[] tmp2 = new int[2];
        tmp1[0] = binominals.get(0).getCoef1();
        tmp1[1] = binominals.get(0).getCoef2();
        tmp2[0] = binominals.get(1).getCoef1();
        tmp2[1] = binominals.get(1).getCoef2();
        int[] accum = polinomMult(tmp1, tmp2);

        for(int i = 2; i < binominals.size(); i++){
            tmp1[0] = binominals.get(i).getCoef1();
            tmp1[1] = binominals.get(i).getCoef2();
            accum = polinomMult(accum, tmp1);
        }
        System.out.println("Accum");
        for(int i = 0; i < accum.length; i++){
            System.out.print(accum[i] + " ");
        }
        System.out.println();
        GeneratingPolynom = new int[D];
        System.out.println("Generic");
        for(int i = 0; i <= D - 1; i++){
            GeneratingPolynom[i] = accum[i];
            System.out.print(GeneratingPolynom[i] + " ");
        }
        System.out.println();
    }

    public int[] codeWorde(){
        int[] inputWord = new int[K];
        int[] codeWordtmp;
        int[] codeWord = new int[N];
        Scanner sc = new Scanner(System.in);
        System.out.println("Please input useful bits:");
        for(int i = 0; i < K; i++){
            inputWord[i] = sc.nextInt();
        }
        System.out.println();
        for(int i = 0; i < inputWord.length; i++){
            inputWord[i]  = mod(inputWord[i]);
        }
        System.out.println("Input word");
        for(int i = 0; i < inputWord.length; i++){
            System.out.print(inputWord[i] + " ");
        }
        int[] tmp1 = new int[N];
        int[] tmp2 = new int[N];
        /*for(int i = 0; i < K; i++){
            tmp1[i] = inputWord[i];
        }
        for(int i = 0; i < tmp1.length; i++){
            tmp2[i] = tmp1[i];
        }*/
        int[] u = new int[N - K + 1];
        u[0] = 1;
        inputWord = polinomMult(inputWord, u);
        for(int i = 0; i < inputWord.length; i++){
            tmp2[i] = inputWord[i];
        }
        int[] ost = polinomDiv(tmp2, GeneratingPolynom);
        System.out.println();
        codeWordtmp = polinomMult(inputWord, GeneratingPolynom);


        for(int i = 0; i < N; i++){
           codeWord[i] = codeWordtmp[i];
        }
        //inputWord = polinomDiv(codeWord, GeneratingPolynom);
        return codeWord;
    }

    public int[] decodeWord(int[] word){
        boolean right = true;
        int[] S = new int[2*t];
        int[] input = new int[word.length];
        for(int i = 0; i < word.length; i++){
            input[i] = word[i];
        }
        int sum = 0;
        //Формирование синдрома
            for (int i = 1; i <= 2*t; i++) {
                sum = 0;
                for (int j = 0; j < word.length; j++) {
                    int alpha = (int) Math.pow(a_p, i);
                    alpha = mod(alpha);
                    int deg = (int) Math.pow(alpha, (word.length - 1) - j);
                    int num = word[j] * deg;
                    S[i-1] = S[i-1] + mod(num);

                }
                S[i-1] = mod(S[i-1]);
                sum = sum + S[i-1];
                if(sum != 0){
                    right = false;
                }
            }
            if(right){
                System.out.println("There is no mistakes");
                return word;
            }
            int[] lock = new int [N];
            int check = 0;
            //Многочлен локаторов ошибок Берликамп-Месси
            /*int[] V = new int[] {0,1,1,1,1,0,1,0,1,1,0,0,1,0,0};
            int[] L = BerlikampMessi(V);
            int[] L = berlecampMessayAlgorithm(V);
            S = new int[]{5, 5, 4, 0};*/
        Polynomial resultPolynomial = berlecampMessayAlgorithm(S, GF);
        int[] result = resultPolynomial.getCoefficients();
        int[] L = new int[result.length];
        int c = result.length - 1;
        for(int i = 0; i < L.length; i++){
            L[c] = result[i];
            c--;
        }
            //Вычисление корней многочлена локаторов ошибок
            for(int i = 0; i < N ; i++){
                int alpha = (int) Math.pow(a_p, i );
                for (int j = 0; j < L.length; j++) {
                    check += L[j] * (int) Math.pow(alpha, (L.length - 1) - j);
                    check = mod(check);
                }
                if(check == 0){
                   lock[i - 1] = 1;
                   /* int alpha_rev = getInverse(alpha, GF);
                    for(int j = 0; j < N; j++){
                        int tmp = (int)Math.pow(a_p, j);
                        tmp = mod(tmp);
                        if(tmp == alpha_rev){
                            lock[j] = 1;
                            break;
                        }
                    }*/
                }
                check = 0;
            }
        System.out.print("L: ");
        for(int i = 0; i < lock.length; i++){
            System.out.print(i + " pos = " + lock[i] + ", ");
        }
        System.out.println();
            /*for(int i = 0; i < lock.length; i++){
                if(lock[i] != 0){
                    int alpha = (int)Math.pow(a_p, i);
                    alpha = mod(alpha);
                    int alpha_rev = getInverse(alpha, GF);
                    lock[i] = 0;
                    for(int j = 0; j < N; j++){
                        int tmp = (int)Math.pow(a_p, j);
                        tmp = mod(tmp);
                        if(tmp == alpha_rev){
                            lock[j] = 1;
                            break;
                        }
                    }
                }
            }*/
            //Вычисление формальной производной
            int[] Lderiv =  polinomDeriv(L);
            //Вычисление многочлена ошибок
            int[] W = polinomMult(S, L);
            int[] x = new int[N - K + 1];
            x[0] = 1;
            W = polinomDiv(W,x);
            int[] e = new int[lock.length];
            int res1 = 0;
            int res2 = 0;
            //Нахождение значений ошибки
        for(int i = 0; i < lock.length; i++){
            if(lock[i] != 0){
                int alpha = (int) Math.pow(a_p, i + 1);
                alpha = mod(alpha);
                int a_rev = getInverse(alpha, GF);
                for(int j = W.length-t; j < W.length; j++){
                    res1 += W[j]*(int)Math.pow(a_rev, (W.length - 1) - j);
                }
                for(int j = 0; j < Lderiv.length; j++){
                    res2 += Lderiv[j]*(int)Math.pow(a_rev, (Lderiv.length - 1) - j);
                }
                res2 = mod(res2);
              /*int newAlpha = getInverse(res2, GF);
              int eps = mod(res2 * res1);*/
                if(trigger <= t) {
                    e[i] = mst[i] * lock[i];
                }
            }
        }
        System.out.print("E: ");
        for(int i = 0; i < e.length; i++){
            System.out.print(i + " pos = " + e[i] + ", ");
        }
        System.out.println();
        for (int i = 0; i < input.length; i++) {
            input[i] = input[i] - e[i];
            input[i] = mod(input[i]);
        }
        Polynomial wrd = new Polynomial(word);
        Polynomial gen = new Polynomial(GeneratingPolynom);
        return input;
    }

    private int findePrimitive(){
        int a = 1;
        boolean isone;
        for (int i = 1; i < GF; i++){
            isone = false;
            int prev = (int) Math.pow(i, (GF  - 1) );
            prev = mod(prev);
            for(int j = 1; j < GF - 1; j++){
                int tmp = (int) Math.pow(i, j);
                tmp = mod(tmp);
                if(tmp == 1){
                    isone = true;
                }
            }
            if(prev == 1 && !isone){
                a = i;
            }
        }
        return a;
    }

    public int[] addMistakes(int[] word){
        int[] f = new int[word.length];
       f[1] = 3;
        f[4] = 2;
       //f[2] = 1;

        mst = f;
        for(int i = 0; i < f.length; i++){
            if(f[i] != 0) {
                trigger += 1;
            }
        }
        for(int i = 0; i < f.length; i++){
            word[i] = mod(word[i] + f[i]);
        }

        return word;
    }
    private int mod (int elem){
        if (elem < GF && elem > 0){
            return elem;
        }
        else{
            while(elem >= GF){
                elem = elem  - GF;
            }
        }
        if(elem < 0){
            while(elem < 0) {
                elem = elem + GF;
            }
        }
        return elem;
    }

    private int[] polinomMult(int[] pol1, int[] pol2){
        int[] resPoly = new int[pol1.length + pol2. length - 1];
        for (int i = pol1.length - 1; i >= 0; i--) {
            for (int j = pol2.length - 1; j >= 0; j--) {
                resPoly[i + j] += pol1[i] * pol2[j];

            }
        }
        for(int i = 0; i < resPoly.length; i++){
            resPoly[i] = mod(resPoly[i]);
        }
        return resPoly;
    }
    private int[] polinomDeriv(int[] pol1){
        int[] deriv = new int[pol1.length];
        for(int i = 0; i < pol1.length; i++){
            int s = pol1.length - i - 1;
            int coef = pol1[i]*s;
            coef = mod(coef);
            int pos = (pol1.length - 1) - (pol1.length - i - 1 - 1);
            if(pos < deriv.length) {
                deriv[pos] = coef;
            }
        }
        int k = 0;
        while(deriv[k] == 0){
            k++;
        }
        int[] res = new int[deriv.length - k];
        for(int i = 0; i < res.length; i++){
            res[i] = deriv[k];
            k++;
        }
        return res;
    }
    private int[] polinomDiv(int[] pol1, int[] pol2){

        int[] qpol;
        int[] respol;
        int k = 0;
        if(pol1.length < pol2.length){
            respol = new int[pol1.length];
            return respol;
        }
        else{
            for(int i = 0; i < pol2.length; i++){
                if(pol2[i] != 0){
                    k = i;
                    break;
                }
            }
        }

        k = 50;

        while(true){
            for (int i = 0; i < pol1.length; i++){
                if(pol1[i] != 0){
                    k = i;
                    break;
                }
            }
            if(k == 50){
                respol = new int[pol1.length];
                return respol;
            }
            if(pol1.length - k  < pol2.length){
                return pol1;
            }
            for(int j = 0; j < pol2.length; j++){
                pol1[k] = pol1[k] - pol2[j];
                if(pol1[k] < 0){
                    while(pol1[k] < 0){
                        pol1[k] += GF;
                    }
                }
                k = k + 1;
            }
        }

    }
    /*private BigInteger algorythmEuclid(BigInteger n, BigInteger p){
        BigInteger Rminus = p;
        BigInteger Rzero = n;
        BigInteger Yminus = BigInteger.ZERO;
        BigInteger Yzero = BigInteger.ONE;
        BigInteger q;
        BigInteger Rnow;
        BigInteger Ynow;
        BigInteger a;
        while(true){
            q = Rminus.divide(Rzero);
            Rnow = Rminus.subtract(Rzero.multiply(q));
            while(Rnow.compareTo(BigInteger.ZERO) < 0){
                Rnow = Rnow.add(p);
            }
            if(Rnow.equals(BigInteger.ZERO)){
                a = Yzero;
                break;
            }
            Ynow = Yminus.subtract(Yzero.multiply(q));
            while(Ynow.compareTo(BigInteger.ZERO) < 0){
                Ynow = Ynow.add(p);
            }
            Rminus = Rzero;
            Rzero = Rnow;
            Yminus = Yzero;
            Yzero = Ynow;

        }
        return a;
    }*/
   /* public int[] BerlikampMessi(int[] S){

        int n = S.length;
        int[] cD = new int[2];
        cD[0] = 0;
        cD[1] = 1;
        int L = 0;
        int m = -1;
        int[] bD = new int[2];
        bD[0] = 0;
        bD[1] = 1;
        int[] tD = new int[1];
        int accum = 0;
        int delta = 0;
        int[] Dpoly;
        int[] C;
        int N_alg = 0;

        while(N_alg < n){
            C = new int[cD.length];
            for(int i = 0; i < cD.length; i++){
                C[i] = cD[i];
            }
            for(int i = 1; i < L + 1; i++){
                accum = accum + S[N_alg - i] * C[i];
            }

            delta = mod(S[N_alg] + accum);
            accum = 0;

            if(delta != 0){
                tD = cD;
                Dpoly = new int[N_alg - m + 1];
                Dpoly[0] = 1;
                int[] tmp = polinomMult(bD, Dpoly);
                for(int i = 0; i < Math.min(cD.length, tmp.length); i++){
                    tmp[i] = mod (tmp[i] + cD[i]);
                }
                cD = tmp;
            }

            if(L <= N_alg/2 && delta != 0){
                L = N_alg + 1 - L;
                m = N_alg;
                bD = tD;
            }
            N_alg++;
        }
        int[] answer = new int[L + 1];
        int k = answer.length - 1;
        for(int i = 0; i < answer.length; i++){
            answer[k] = cD[i];
            k--;
        }
        return answer;
    }*/

    private static Polynomial berlecampMessayAlgorithm(int[] s, int p) {
        Polynomial c = new Polynomial(new int[]{1});  // result register properties
        Polynomial b = new Polynomial(new int[]{1});

        int L = 0;
        int z = 1;
        int n = 0;
        int delta = 1;

        while(n < s.length) {
            int d = s[n];
            for (int i = 1; i <= L; i++) {
                d += c.getCoefficients()[i] * s[n - i];
            }
            d %= p;

            if (d != 0) {
                Polynomial t = c;

                int deltaInverse = getInverse(delta, p);
                Polynomial tmp = Polynomial.getOneVariablePolynomial(d * deltaInverse, z);
                Polynomial mult = Polynomial.multiply(b,tmp, p);
                c = new Polynomial(Polynomial.subtract(c, mult, p).getCoefficients());

                if (2 * L <= n) {
                    L = n + 1 - L;
                    z = 1;
                    b = t;
                    delta = d;
                } else {
                    z++;
                }
            } else {
                z++;
            }

            n++;
        }

        return c;
    }

    public static int getInverse(int a, int m) {
        int m0 = m;
        int y = 1, x = 0;

        while (m > 1)
        {
            int q = m / a;
            int t = a;

            a = m % a;
            m = t;
            t = y;

            y = x - q * y;
            x = t;
        }

        if (x < 0)
            x += m0;

        return x;
    }



    public int[] getGeneratingPolynom(){
        return GeneratingPolynom;
    }

    public void showGenericPolynom(){
        for(int i = 0; i < GeneratingPolynom.length; i++){
            int deg = GeneratingPolynom.length - i - 1;
            System.out.print(GeneratingPolynom[i]+"x^" + deg);
            if(i < GeneratingPolynom.length - 1){
                System.out.print(" + ");
            }
        }
        System.out.println();
    }




























}

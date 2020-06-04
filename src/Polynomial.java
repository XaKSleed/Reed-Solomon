public class Polynomial {
    private int[] coefficients;

    public Polynomial(int[] coefficients) {
        this.coefficients = coefficients;
    }

    public static Polynomial multiply(Polynomial p1, Polynomial p2, int modulo) {
        int[] resultCoefficients = new int[p1.length() + p2.length() - 1];

        for(int i = 0; i < p1.length(); i++) {
            for(int j = 0; j < p2.length(); j++) {
                resultCoefficients[i + j] += p1.getCoefficients()[i] * p2.getCoefficients()[j];
                resultCoefficients[i + j] %= modulo;
            }
        }

        return new Polynomial(resultCoefficients);
    }

    public static Polynomial subtract(Polynomial p1, Polynomial p2, int modulo) {
        int len = Math.max(p1.length(), p2.length());
        int[] resultCoefficients = new int[len];

        for(int i = 0; i < len; i++) {
            if(i >= p1.length()) {
                int sub = (-p2.getCoefficients()[i]) + modulo;
                resultCoefficients[i] = sub % modulo;
            } else if(i >= p2.length()) {
                resultCoefficients[i] = p1.getCoefficients()[i];
            } else {
                int sub = (p1.getCoefficients()[i] - p2.getCoefficients()[i]) + modulo;
                resultCoefficients[i] = sub % modulo;
            }
        }

        return new Polynomial(resultCoefficients);
    }

    public static Polynomial getOneVariablePolynomial(int coefficient, int deg) {
        int[] coefs = new int[deg + 1];
        coefs[deg] = coefficient;
        return new Polynomial(coefs);
    }

    public int length() {
        return this.coefficients.length;
    }

    public int[] getCoefficients() {
        return coefficients;
    }

    public void multiplyNumber(int n, int modulo) {
        for(int i = 0; i < this.coefficients.length; i++) {
            this.coefficients[i] = (this.coefficients[i] * n) % modulo;
        }
    }

    @Override
    public String toString() {
        if(this.coefficients.length == 0) { return ""; }

        StringBuilder result = new StringBuilder("");

        for(int i = 0; i < this.coefficients.length; i++) {
            if(this.coefficients[i] != 0) {
                if(!result.toString().equals("")) {
                    result.append(" + ");
                }

                if(i == 0) {
                    result.append(this.coefficients[i]);
                    continue;
                }

                if (this.coefficients[i] == 1) {
                    result.append("x");
                } else {
                    result.append(this.coefficients[i]).append("x");
                }

                if (i > 1) {
                    result.append("^").append(i);
                }
            }
        }

        return result.toString();
    }
}
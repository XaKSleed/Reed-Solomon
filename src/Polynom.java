
import java.util.Arrays;

/**
 * Created by Elena on 16.02.2017
 */

public class Polynom {
    private double array[];
    private int deg;

    public int gettdeg() {
        return array.length-1;
    }

    public double getCoefficient(int deg) {
        if (deg < 0 || deg > this.gettdeg())
            throw new IllegalArgumentException("Неверная степень");
        return this.array[deg];
    }

    @Override
    public String toString() {
        StringBuilder s = new StringBuilder("");
        for (int i=this.deg; i>=0; i--) {
            if (this.array[i] == 0) continue;
            if ((this.array[i] > 0) && (i!=this.deg))
                s.append("+");
            if ((int)this.array[i] == this.array[i]) s.append((int)this.array[i]);
            else s.append(this.array[i]);
            if (i == 1) s.append("x");
            else if (i != 0) s.append("x^" + i);
        }
        return s.toString();
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        Polynom polynom = (Polynom) o;

        if (deg != polynom.deg) return false;
        return Arrays.equals(array, polynom.array);

    }

    @Override
    public int hashCode() {
        int result = deg;
        result = 31 * result + Arrays.hashCode(array);
        return result;
    }

    public Polynom(double[] array) {
        this.array = new double [array.length];
        for (int i=0;i<this.array.length;i++) {
            this.array[i] = array[array.length-1-i];
        }
        this.deg = gettdeg();
    }

    private Polynom(int deg) {
        this.deg = deg;
        this.array = new double[this.deg+1];
    }

    public Polynom(Polynom p) {
        this.deg = p.deg;
        this.array = new double[this.deg+1];
        System.arraycopy(p.array, 0, this.array, 0, p.deg+1);
    }

    public void show() {
        System.out.println(this);
    }


    public Polynom sum(Polynom p2) {
        int max = Math.max(this.deg, p2.deg), min = Math.min(this.deg, p2.deg);
        Polynom res = new Polynom(max); //новый полином размера большего из складываемых полиномов
        for (int i=0; i<=min; i++) {
            res.array[i] = this.array[i] + p2.array[i];
        }
        if (this.deg > p2.deg) System.arraycopy(this.array, min+1, res.array, min+1, max-min);
        else if (p2.deg > this.deg) System.arraycopy(p2.array, min+1, res.array, min+1, max-min);
        res.cut0();
        return res;
    }

    public Polynom invert() {
        Polynom res = new Polynom(this);
        for (int i=0; i<=res.deg; i++)
            res.array[i] = -res.array[i];
        return res;
    }

    public Polynom dif(Polynom p2) {
        Polynom copyp2 = p2.invert();
        Polynom res = this.sum(copyp2);
        return res;
    }

    public double whatx(int x) {
        double result = this.array[0];
        for (int i=1;i<=this.deg;i++) {
            result += this.array[i]*Math.pow((double)x, (double)i);
        }
        return result;
    }

    private void cut0() {
        if (this.deg == 0) return;
        int i;
        for (i = this.deg; i >= 0; i--) {
            if (this.array[i] != 0) break;
        }
        if (i == this.deg || i == -1) return;
        double newar[] = new double[this.deg+1];
        System.arraycopy(this.array, 0, newar, 0, this.deg+1);
        this.deg = i;
        this.array = new double[this.deg + 1];
        System.arraycopy(newar, 0, this.array, 0, this.deg + 1);
        return;
    }

    public Polynom multipl(Polynom p2){
        Polynom res = new Polynom(this.deg+p2.deg);
        for (int i=this.deg;i>=0;i--) {
            for (int j=p2.deg;j>=0;j--) {
                res.array[i+j] += this.array[i] * p2.array[j];
            }
        }
        res.cut0();
        return res;
    }

    public Polynom div (Polynom p2) {
        if (this.equals(p2)) {
            Polynom res = new Polynom(1);
            res.array[0] = 1;
            return res;
        }
        this.cut0(); p2.cut0();
        if (p2.array[p2.deg] == 0.0) throw new IllegalArgumentException("На ноль делить нельзя!");
        Polynom ost = new Polynom(this);
        Polynom res = new Polynom(this.deg);
        Polynom multipleres = new Polynom(p2.deg);
        Polynom mnog = new Polynom(this.deg);
        double mn;
        while (ost.deg >= p2.deg) {
            mn = ost.array[ost.deg] / p2.array[p2.deg];
            res.array[ost.deg - p2.deg] = mn;
            Arrays.fill(mnog.array, 0);
            mnog.array[ost.deg - p2.deg] = mn;
            multipleres = p2.multipl(mnog);
            ost = ost.dif(multipleres);
        }
        res.cut0();
        return res;
    }

    public Polynom divost (Polynom p2) {
        if (this.equals(p2)) {
            Polynom ost = new Polynom(1);
            ost.array[0] = 0;
            return ost;
        }
        this.cut0(); p2.cut0();
        if (p2.array[p2.deg] == 0.0) throw new IllegalArgumentException("На ноль делить нельзя!");
        Polynom ost = new Polynom(this);
        Polynom multipleres = new Polynom(p2.deg);
        Polynom mnog = new Polynom(this.deg);
        double mn;
        while (ost.deg >= p2.deg) {
            mn = ost.array[ost.deg] / p2.array[p2.deg];
            Arrays.fill(mnog.array, 0);
            mnog.array[ost.deg - p2.deg] = mn;
            multipleres = p2.multipl(mnog);
            ost = ost.dif(multipleres);
        }
        ost.cut0();
        return ost;
    }
}
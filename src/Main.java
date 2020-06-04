public class Main {

    public static void main(String[] args) {

        ReedSolomonCode code1 = new ReedSolomonCode(6, 5);
        code1.showGenericPolynom();
       int[] codeword1 = code1.codeWorde();
        System.out.println("Code word");
       for(int i = 0; i < codeword1.length; i++){
           System.out.print(codeword1[i] + " ");
        }
        System.out.println();
        code1.getGeneratingPolynom();
   int[] mistword = code1.addMistakes(codeword1);
        System.out.println("Mistaken");
        for(int i = 0; i < mistword.length; i++){
            System.out.print(mistword[i] + " ");
        }
        System.out.println();

        System.out.println("Decode");
        int[] res = code1.decodeWord(codeword1);
        System.out.println("Correct word");
        for(int i = 0; i < res.length; i++){
            System.out.print(res[i] + " ");
        }
        System.out.println();


    }
}

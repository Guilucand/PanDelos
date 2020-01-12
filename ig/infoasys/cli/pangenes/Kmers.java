package infoasys.cli.pangenes;

public class Kmers {

    static long defaultStep = ('Z' - 'A' + 1);

    static long getProteinMapping(char base) {
        switch (base) {
            case 'A': return 0;
            case 'C': return 1;
            case 'G': return 2;
            case 'T': return 3;
            case 'N': return 4;
        }
        return -1;
    }

    static long getDefaultMapping(char base) {
        return base - 'A';
    }




    static long computeNext(long current, char value, long kpow) {
        current *= defaultStep;
        current += getDefaultMapping(value);
        return current % kpow;
    }

    public static long[] computeFromSequence(String seq, int k) {
        long kpow = 1;
        for(int i = 0; i < k; i++) kpow *= defaultStep;

        long[] results = new long[seq.length() - k + 1];

        long rank = 0;
        for (int i = 0; i < k; i++) {
            rank = computeNext(rank, seq.charAt(i), kpow);
        }
        results[0] = rank;

        for (int i = k; i < seq.length(); i++) {
            rank = computeNext(rank, seq.charAt(i), kpow);
            results[i-k+1] = rank;
        }
        return results;
    }

}

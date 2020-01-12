package infoasys.cli.pangenes;

public class PangeneNative {

    PangeneNative(int k, PangeneIData data) {
        initialize(k);
        preprocessSequences(data);
    }

    private native void initialize(int kvalue);

    private native void preprocessSequences(PangeneIData data);

    private native long computePair(int[] genome1, int[] genome2);

    private native void computeSequenceScores(long ref, int row_start, int row_end, PairScores out_scores);

    private native void freePairStruct(long structure);

    public PairInfo generatePairInfo(int[] genome1, int[] genome2) {
        return new PairInfo(genome1, genome2);
    }

    public class PairInfo {
        private long native_addr;
        private int len;

        protected PairInfo(int[] genome1, int[] genome2) {
            native_addr = computePair(genome1, genome2);
            len = genome1.length + genome2.length;
        }

        public PairScores computePairScoresAndFree() {
//            System.out.println("CMPTT " + native_addr);
            PairScores result = new PairScores();
            computeSequenceScores(native_addr, 0, len, result);
            freePairStruct(native_addr);
//            System.out.println("CMPTD " + native_addr);
            return result;
        }
    }
}

package infoasys.cli.pangenes;

public class PangeneNative {

    PangeneNative(int k, PangeneIData data) {
        preprocessSequences(data, k, false);
    }

    private PangeneNative() {}
    static void printComplexity(int k, PangeneIData data) {
        new PangeneNative().preprocessSequences(data, k, true);
    }

    private native void preprocessSequences(PangeneIData data, int kvalue, boolean onlyComplexity);
    private native void computeScores(int genome, Scores out_scores, int step_size);

    public Scores generateScoresPart(int genome, boolean multithread) {
        Scores result = new Scores();
        computeScores(genome, result, multithread ? 2048 : Integer.MAX_VALUE);
        return result;
    }
}

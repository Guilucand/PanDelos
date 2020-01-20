package infoasys.cli.pangenes;

public class PangeneNative {

    PangeneNative(int k, PangeneIData data) {
        initialize(k);
        preprocessSequences(data);
    }

    private native void initialize(int kvalue);
    private native void preprocessSequences(PangeneIData data);
    private native void computeScores(int genome, Scores out_scores);

    public Scores generateScoresPart(int genome) {
        Scores result = new Scores();
        computeScores(genome, result);
        return result;
    }
}

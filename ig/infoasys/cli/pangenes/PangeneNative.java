package infoasys.cli.pangenes;

public class PangeneNative {

    PangeneNative(int k, PangeneIData data) {
        preprocessSequences(data, k);
    }

    private native void preprocessSequences(PangeneIData data, int kvalue);
    private native void computeScores(int genome, Scores out_scores);

    public Scores generateScoresPart(int genome) {
        Scores result = new Scores();
        computeScores(genome, result);
        return result;
    }
}

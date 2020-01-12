package infoasys.cli.pangenes;

public class PairScores {
    int scoresCount;
    // For each element
    float[] scores;
    // For each element
    float[] percs;
    // For each element
    float[] tr_percs;
    // For each element
    int[] row;
    // For each element
    int[] column;
    // For each element
    int[] firstSeqIndex;
    // For each element
    int[] secondSeqIndex;
    // For each element
    boolean[] sameGenome;

    // For each row
    float[] max_intra_score;
    // For each row
    float[] max_inter_score;
}


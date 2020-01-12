package infoasys.cli.pangenes;

import java.io.FileWriter;
import java.io.IOException;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.Vector;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;

import infoasys.core.common.util.Pair;

public class Pangenes {

	public static void usage() {
		System.err.println("Usage: cmd ifile.faa k ofile.net");
	}

	public static void main(String[] args) {

		System.out.println("Working Directory = " +
				System.getProperty("user.dir"));

		System.loadLibrary("native");

		String ifile = null;
		int k = 1;
		String netofile = null;


		try {
			ifile = args[0];
			k = Integer.parseInt(args[1]);
			netofile = args[2];
		} catch (Exception e) {
			usage();
			System.exit(1);
		}

		float kk = (float) k;

		PangeneIData pid;
		try {
			pid = PangeneIData.readFromFile(ifile);
		} catch (Exception e) {
			e.printStackTrace();
			return;
		}

		long startt = System.currentTimeMillis();
		PangeneNative nativ = new PangeneNative(k, pid);


		int nofGenomes = pid.genomeNames.size();
		Map<Integer, Vector<Integer>> genomeSets = pid.getGenomeSets();
		PangeneNet pnet = new PangeneNet();

		/*
		 * Total complexity: |gset|**2 * ((|seq(gi)| + |seq(gj)|)**2 + ...)
		 * Worst case complexity: |gset|**2 * ((|seq(gi)| + |seq(gj)|)**2 + len(seq(gi)) * (|seq(gi)| + |seq(gj)|))
		 * New complexity: |gset|**2 * (|seq(gi)| * |seq(gj)| * len(seq(gi)))
		 *
		 * */
		ExecutorService threadService = Executors.newFixedThreadPool(12);

		for (int g1 = 0; g1 < nofGenomes; g1++) {
			for (int g2 = g1 + 1; g2 < nofGenomes; g2++) {

				final int g1f = g1;
				final int g2f = g2;

				threadService.submit(() -> {
					System.out.println("##\t" + g1f + "\t" + g2f);

					int[] allowedFirst = genomeSets.get(g1f).stream().mapToInt(i -> i).toArray();
					int[] allowedSecond = genomeSets.get(g2f).stream().mapToInt(i -> i).toArray();

					int sequencesCount = allowedFirst.length + allowedSecond.length;

					PangeneNative.PairInfo pairInfo = nativ.generatePairInfo(allowedFirst, allowedSecond);
					PairScores pairScores = pairInfo.computePairScoresAndFree();

					if (g1f == g2f) {
						/* connect identical sequences of same genome */
						for (int i = 0; i < pairScores.scoresCount; i++) {
							if ((pairScores.row[i] < pairScores.column[i]) && pairScores.scores[i] == 1.0) {
								pnet.addConnection(pairScores.firstSeqIndex[i], pairScores.secondSeqIndex[i], pairScores.scores[i]);
							}
						}
					}
					else {
						boolean[] engagged = new boolean[sequencesCount];
						Arrays.fill(engagged, false);

						/* get inter bbh */
						/* also, calcolate threshold for inter non-bbh as the average non-null and non-bbh scores*/
						float inter_thr_sum = 0.0f;
						float inter_thr_count = 0.0f;
						float inter_max_score = 0.0f;
						float inter_min_score = 1.0f;

						float inter_max_perc = 0.0f;
						float inter_min_perc = 1.0f;

						FileWriter debugtest = null;

						try {
							debugtest = new FileWriter("dtest.txt");
						} catch (IOException e) {
							e.printStackTrace();
						}

						//0-
						//x0
						for (int i = 0; i < pairScores.scoresCount; i++) {
							int row = pairScores.row[i];
							int column = pairScores.column[i];
							float score = pairScores.scores[i];
							int firstSeqIndex = pairScores.firstSeqIndex[i];
							int secondSeqIndex = pairScores.secondSeqIndex[i];

							float perc = pairScores.percs[i];
							float otherPerc = pairScores.tr_percs[i];

							if (!pairScores.sameGenome[i] && pairScores.max_inter_score[row] > 0.0f && pairScores.max_inter_score[column] > 0.0f) {

//								OTHER SCORE IS ALWAYS EQUAL TO CURRENT SCORE!!
//								float otherScore = (otherIndex >= 0) ? pairScores.scores[otherIndex] : 0.0f;
								float otherScore = score;

								if ((score == pairScores.max_inter_score[row] && otherScore == pairScores.max_inter_score[column])) {
									pnet.addConnection(firstSeqIndex, secondSeqIndex, score);
									pnet.addConnection(secondSeqIndex, firstSeqIndex, otherScore);
									engagged[row] = true;
									engagged[column] = true;

									inter_thr_sum += score + otherScore;
									inter_thr_count += 2;
									if (score < 1.0 && score > inter_max_score) {
										inter_max_score = score;
									} else if (otherScore < 1.0 && otherScore > inter_max_score) {
										inter_max_score = otherScore;
									}

									if (score > 0.0 && score < inter_min_score) {
										inter_min_score = score;
									} else if (otherScore > 0.0 && otherScore < inter_min_score) {
										inter_min_score = otherScore;
									}

									inter_max_perc = Math.max(inter_max_perc, Math.max(perc, otherPerc));
									inter_min_perc = Math.min(inter_min_perc, Math.min(perc, otherPerc));
								}
							}
						}

						/* calcolate threshold for inter non-bbh */
						float inter_thr = inter_thr_sum / inter_thr_count;

						System.out.println("score\t" + inter_thr + "\t" + inter_min_score + "\t" + inter_max_score);
						System.out.println("perc\t" + inter_min_perc + "\t" + inter_max_perc);
						System.out.println(pnet.countNodes() + "\t" + pnet.countEdges());

						//-x
						//x0
						for (int i = 0; i < pairScores.scoresCount; i++) {
							if (pairScores.sameGenome[i]) {
								float score = pairScores.scores[i];
								int row = pairScores.row[i];
								int column = pairScores.column[i];
								int firstSeq = pairScores.firstSeqIndex[i];
								int secondSeq = pairScores.secondSeqIndex[i];
								if (engagged[row]) {
									/* get identical paralogs */
									if (score == 1.0) {
										pnet.addConnection(firstSeq, secondSeq, score);
									}
									/* get paralog / intra bbh, filtered by inter max scores */
									else if (score == pairScores.max_intra_score[row] &&
											score == pairScores.max_intra_score[column] &&
											score >= inter_max_score) {
										try {
											debugtest.write(firstSeq + " " + secondSeq + " " + score + "\n");
										} catch (IOException e) {
											e.printStackTrace();
										}
										pnet.addConnection(firstSeq, secondSeq, score);
									}
								}
							}
						}
						System.out.println(pnet.countNodes() + "\t" + pnet.countEdges());
					}
				});
			}
		}

		boolean terminated = false;
		while (!terminated) {
			try {
				threadService.shutdown();
				threadService.awaitTermination(1, TimeUnit.SECONDS);
				terminated = threadService.isTerminated();
			}
			catch (Exception ignored) {
			}
		}
		long tott = System.currentTimeMillis() - startt;
		System.out.println("Total time: " + tott);


		System.out.println("----------");
		System.out.println("undirected degree distribution");
		pnet.printUndirectedDegreeDistribution();

		System.out.println("----------");
		System.out.println("directed degree distribution");
		pnet.printDegreeDistribution();

		System.out.println("----------");
		Map<Integer, Integer> cocodistr = new TreeMap<Integer, Integer>();
		List<Set<Integer>> cocos = pnet.undirectedConnectedComponets();
		for (Set<Integer> coco : cocos) {
			cocodistr.put(coco.size(), cocodistr.getOrDefault(coco.size(), 0) + 1);
		}
		System.out.println("----------");
		System.out.println("CoCo sizes");
		for (Map.Entry<Integer, Integer> en : cocodistr.entrySet()) {
			System.out.println(en.getKey() + "\t" + en.getValue());
		}


		System.out.println("----------");
		System.out.println("writing into " + netofile + "");
		try {
			pnet.saveToFile(netofile, false);
		} catch (Exception e) {
			e.printStackTrace();
			return;
		}
	}
}

package infoasys.cli.pangenes;

import java.util.Arrays;
import java.util.BitSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.atomic.AtomicInteger;

public class Pangenes {

	public static void main(String[] args) {

		System.out.println("Working Directory = " +
				System.getProperty("user.dir"));

		System.loadLibrary("native");

		Cli cli = new Cli(args);

		PangeneIData pid;
		try {
			pid = PangeneIData.readFromFile(cli.InputFile);
		} catch (Exception e) {
			e.printStackTrace();
			return;
		}

		if (cli.OnlyOps) {
			PangeneNative.printComplexity(cli.KValue, pid);
			return;
		}

		long startt = System.currentTimeMillis();
		PangeneNative nativ = new PangeneNative(cli.KValue, pid);

		AtomicInteger doneGenomes = new AtomicInteger();

		int nofGenomes = pid.genomeNames.size();
		PangeneNet pnet = new PangeneNet();

		/*
		 * Total complexity: |gset|**2 * ((|seq(gi)| + |seq(gj)|)**2 + ...)
		 * Worst case complexity: |gset|**2 * ((|seq(gi)| + |seq(gj)|)**2 + len(seq(gi)) * (|seq(gi)| + |seq(gj)|))
		 * New complexity: |gset|**2 * (|seq(gi)| * |seq(gj)| * len(avg(seq(gi), seq(g2))))
		 * New complexity2: (|seq| * |seq|) + |seq| * avg(kmers_repeats)))
		 *
		 * */

		ExecutorService threadService = Executors.newFixedThreadPool(cli.ThreadsNum);

		int sequencesCount = pid.sequences.size();

		Object printLock = new Object();

		for (int g = 0; g < nofGenomes; g++) {

			final int finalG = g;

			threadService.submit(() -> {
				System.out.println("Working on genome " + finalG + "/" + nofGenomes);
				final Scores scoresPart = nativ.generateScoresPart(finalG, cli.ThreadsNum != 1);
				System.out.println("Preprocessed genome " + finalG + "/" + nofGenomes);
				System.out.println("Filtered count: " + scoresPart.scoresCount);

				final float[] inter_thr_sum = new float[nofGenomes];
				final float[] inter_thr_count = new float[nofGenomes];
				final float[] inter_max_score = new float[nofGenomes];
				final float[] inter_min_score = new float[nofGenomes];
				Arrays.fill(inter_min_score, 1.0f);

				final float[] inter_max_perc = new float[nofGenomes];
				final float[] inter_min_perc = new float[nofGenomes];
				Arrays.fill(inter_min_perc, 1.0f);

				// The engagged array should be different for each comparison genome
				// leading to an array of [sequencesCount][nofGenomes] size.
				// this should then be used to determine if we should consider
				// the inter_min_score of the i-th genome to add this sequence to the map
				// An alternative way is to iterate two times the scores, the first time updating only
				// the inter_min/max_score arrays, then updating the relative threshold only when we set engagged to true
//				final boolean[] engagged = new boolean[sequencesCount];
				final BitSet[] engagged = new BitSet[sequencesCount];
				for (int i = 0; i < engagged.length; i++) {
					engagged[i] = new BitSet(nofGenomes);
				}

				Arrays.fill(inter_max_perc, 1.0f);

				float min_inter_max_score = 1.0f;

				boolean[] shouldAddConnection = new boolean[scoresPart.scoresCount];

				for (int i = 0; i < scoresPart.scoresCount; i++) {
					if (scoresPart.first_seq_genome[i] != scoresPart.second_seq_genome[i]) {
						if (scoresPart.scores[i] == scoresPart.max_genome_score[scoresPart.scoresMaxMappings[scoresPart.row[i]]][scoresPart.second_seq_genome[i]] &&
								scoresPart.scores[i] == scoresPart.max_genome_score_col[scoresPart.column[i]]) {

							pnet.addConnection(scoresPart.row[i], scoresPart.column[i], scoresPart.scores[i]);
							pnet.addConnection(scoresPart.column[i], scoresPart.row[i], scoresPart.scores[i]);

							shouldAddConnection[i] = true;

							int sg = scoresPart.second_seq_genome[i];
							float score = scoresPart.scores[i];
							float perc = scoresPart.percs[i];
							float otherPerc = scoresPart.tr_percs[i];

							engagged[scoresPart.row[i]].set(sg);
							inter_thr_sum[sg] += 2 * scoresPart.scores[i];
							inter_thr_count[sg] += 2;
							if (score < 1.0 && score > inter_max_score[sg]) {
								inter_max_score[sg] = score;
							}

							if (score > 0.0 && score < inter_min_score[sg]) {
								inter_min_score[sg] = score;
							}

							inter_max_perc[sg] = Math.max(inter_max_perc[sg], Math.max(perc, otherPerc));
							inter_min_perc[sg] = Math.min(inter_min_perc[sg], Math.min(perc, otherPerc));
						}
					}
				}

				for (int i = 0; i < nofGenomes; i++) {

					if (i == finalG) continue;

					float inter_thr = inter_thr_sum[i] / inter_thr_count[i];

					synchronized (printLock) {
						System.out.println("Comparing genome " + finalG + " with " + i + ":");

						System.out.println("Score\t" + inter_thr + "\t" + inter_min_score[i] + "\t" + inter_max_score[i]);
						System.out.println("Perc\t" + inter_min_perc[i] + "\t" + inter_max_perc[i]);
						System.out.println(pnet.countNodes() + "\t" + pnet.countEdges());
					}
				}


				float[] scoresRowThreshold = new float[sequencesCount];
				Arrays.fill(scoresRowThreshold, Float.POSITIVE_INFINITY);

				for (int i = 0; i < scoresPart.scoresCount; i++) {
					if (shouldAddConnection[i]) {
						int row = scoresPart.row[i];
						int sg = scoresPart.second_seq_genome[i];
						scoresRowThreshold[row] = Math.min(scoresRowThreshold[row], inter_max_score[sg]);
					}
				}

//				for (int i = 0; i < nofGenomes; i++) {
//					if (i == finalG) continue; // Exclude current genome from computation
//					min_inter_max_score = Math.min(min_inter_max_score, inter_max_score[i]);
//				}

				/* get inter bbh */
				/* also, calcolate threshold for inter non-bbh as the average non-null and non-bbh scores*/
				for (int i = 0; i < scoresPart.scoresCount; i++) {
//					if (scoresPart.row[i] == 1915 && scoresPart.column[i] == 1921) {
//						int x = 5;
//					}

					 if ((scoresPart.row[i] < scoresPart.column[i]) &&
//							 engagged[scoresPart.row[i]].get(scoresPart.second_seq_genome[i]) &&
							 scoresPart.first_seq_genome[i] == scoresPart.second_seq_genome[i] &&
							 (scoresPart.scores[i] == scoresPart.max_genome_score[scoresPart.scoresMaxMappings[scoresPart.row[i]]][scoresPart.second_seq_genome[i]] &&
									scoresPart.scores[i] == scoresPart.max_genome_score[scoresPart.scoresMaxMappings[scoresPart.column[i]]][scoresPart.second_seq_genome[i]] &&
									scoresPart.scores[i] >= scoresRowThreshold[scoresPart.row[i]]
//									scoresPart.scores[i] >= min_inter_max_score
							)) {

						pnet.addConnection(scoresPart.row[i], scoresPart.column[i], scoresPart.scores[i]);
					}
				}
				doneGenomes.getAndIncrement();

				synchronized (printLock) {
					System.out.println("Done genome " + doneGenomes.get() + "/" + nofGenomes + " => " + scoresPart.scoresCount);
				}
			});
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
		long memUsed = Runtime.getRuntime().totalMemory();
		System.out.println("JVM Memory used: " + (double)memUsed / (1024 * 1024) + "MB");
		System.out.println("Total time: " + tott);

		System.out.println("----------");
		System.out.println("undirected degree distribution");
		pnet.printUndirectedDegreeDistribution();

		System.out.println("----------");
		System.out.println("directed degree distribution");
		pnet.printDegreeDistribution();

		System.out.println("----------");
		Map<Integer, Integer> cocodistr = new TreeMap<>();
		List<Set<Integer>> cocos = pnet.undirectedConnectedComponents();
		for (Set<Integer> coco : cocos) {
			cocodistr.put(coco.size(), cocodistr.getOrDefault(coco.size(), 0) + 1);
		}
		System.out.println("----------");
		System.out.println("CoCo sizes");
		for (Map.Entry<Integer, Integer> en : cocodistr.entrySet()) {
			System.out.println(en.getKey() + "\t" + en.getValue());
		}

		System.out.println("----------");
		System.out.println("writing into " + cli.OutputNetFile + "");
		try {
			pnet.saveToFile(cli.OutputNetFile, false);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}

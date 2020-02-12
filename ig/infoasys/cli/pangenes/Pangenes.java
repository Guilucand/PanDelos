package infoasys.cli.pangenes;

import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.atomic.AtomicInteger;

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

		int limit_rows = 0;

		try {
			ifile = args[0];
			k = Integer.parseInt(args[1]);
			netofile = args[2];
			try {
				limit_rows = Integer.parseInt(args[3]);
			}
			catch (Exception e) {
				limit_rows = Integer.MAX_VALUE;
			}
		} catch (Exception e) {
			usage();
			System.exit(1);
		}

		PangeneIData pid;
		try {
			pid = PangeneIData.readFromFile(ifile, limit_rows);
		} catch (Exception e) {
			e.printStackTrace();
			return;
		}

		AtomicInteger doneGenomes = new AtomicInteger();

		long startt = System.currentTimeMillis();
		PangeneNative nativ = new PangeneNative(k, pid);


		int nofGenomes = pid.genomeNames.size();
		PangeneNet pnet = new PangeneNet();

		/*
		 * Total complexity: |gset|**2 * ((|seq(gi)| + |seq(gj)|)**2 + ...)
		 * Worst case complexity: |gset|**2 * ((|seq(gi)| + |seq(gj)|)**2 + len(seq(gi)) * (|seq(gi)| + |seq(gj)|))
		 * New complexity: |gset|**2 * (|seq(gi)| * |seq(gj)| * len(avg(seq(gi), seq(g2))))
		 * New complexity2: (|seq| * |seq|) + |seq| * avg(kmers_repeats)))
		 *
		 * */

		int cores = 5;//Runtime.getRuntime().availableProcessors();
		ExecutorService threadService = Executors.newFixedThreadPool(cores);

		int sequencesCount = pid.sequences.size();

		long memUsed = Runtime.getRuntime().totalMemory();
		System.out.println("Memory used: " + (double)memUsed / (1024 * 1024) + "MB");

		Object printLock = new Object();

		for (int g = 0; g < nofGenomes; g++) {

			final int finalG = g;

			threadService.submit(() -> {
				System.out.println("Working on genome " + finalG + "/" + nofGenomes);
				final Scores scoresPart = nativ.generateScoresPart(finalG);
				System.out.println("Preprocessed genome " + finalG + "/" + nofGenomes);

				final float[] inter_thr_sum = new float[nofGenomes];
				final float[] inter_thr_count = new float[nofGenomes];
				final float[] inter_max_score = new float[nofGenomes];
				final float[] inter_min_score = new float[nofGenomes];
				Arrays.fill(inter_min_score, 1.0f);

				final float[] inter_max_perc = new float[nofGenomes];
				final float[] inter_min_perc = new float[nofGenomes];
				Arrays.fill(inter_min_perc, 1.0f);
				final boolean[] engagged = new boolean[sequencesCount];

				Arrays.fill(inter_max_perc, 1.0f);

				float min_inter_max_score = 1.0f;

				for (int i = 0; i < scoresPart.scoresCount; i++) {
					if (scoresPart.first_seq_genome[i] != scoresPart.second_seq_genome[i]) {
						if (scoresPart.scores[i] == scoresPart.max_genome_score[scoresPart.scoresMaxMappings[scoresPart.row[i]]][scoresPart.second_seq_genome[i]] &&
								scoresPart.scores[i] == scoresPart.max_genome_score_col[scoresPart.column[i]]) {
							pnet.addConnection(scoresPart.row[i], scoresPart.column[i], scoresPart.scores[i]);
							pnet.addConnection(scoresPart.column[i], scoresPart.row[i], scoresPart.scores[i]);
							engagged[scoresPart.row[i]] = true;

							int sg = scoresPart.second_seq_genome[i];
							float score = scoresPart.scores[i];
							float perc = scoresPart.percs[i];
							float otherPerc = scoresPart.tr_percs[i];

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

						System.out.println("score\t" + inter_thr + "\t" + inter_min_score[i] + "\t" + inter_max_score[i]);
						System.out.println("perc\t" + inter_min_perc[i] + "\t" + inter_max_perc[i]);
						System.out.println(pnet.countNodes());// + "\t" + pnet.countEdges());
					}
				}

				for (int i = 0; i < nofGenomes; i++) {
					if (i == finalG) continue; // Exclude current genome from computation
					min_inter_max_score = Math.min(min_inter_max_score, inter_max_score[i]);
				}

				/* get inter bbh */
				/* also, calcolate threshold for inter non-bbh as the average non-null and non-bbh scores*/
				for (int i = 0; i < scoresPart.scoresCount; i++) {
					 if ((scoresPart.row[i] < scoresPart.column[i]) &&
							 engagged[scoresPart.row[i]] && scoresPart.first_seq_genome[i] == scoresPart.second_seq_genome[i] &&
							 (scoresPart.scores[i] == 1.0f ||
							(scoresPart.scores[i] == scoresPart.max_genome_score[scoresPart.scoresMaxMappings[scoresPart.row[i]]][scoresPart.second_seq_genome[i]] &&
									scoresPart.scores[i] == scoresPart.max_genome_score[scoresPart.scoresMaxMappings[scoresPart.column[i]]][scoresPart.second_seq_genome[i]] &&
									scoresPart.scores[i] >= min_inter_max_score))) {
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
		System.out.println("Total time: " + tott);
		System.out.println("FMT;" + limit_rows + ";" + tott);


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
		System.out.println("writing into " + netofile + "");
		try {
			pnet.saveToFile(netofile, false);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}

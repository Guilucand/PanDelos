package infoasys.cli.pangenes;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;

public class PangeneNet {

	public static class Edge implements Comparable<Edge> {
		public int node;
		public double score;

		public Edge(int node, double score) {
			this.node = node;
			this.score = score;
		}

		@Override
		public int compareTo(Edge o) {
			return this.node - o.node;
		}

		@Override
		public int hashCode() {
			return this.node;
		}

	}

	public final Map<Integer, Set<Edge>> adjacents;
	Set<Integer> nodes;

	long totalEdges;
	long totalNodes;

	public PangeneNet() {
		this.adjacents = new HashMap<>();
		this.nodes = new HashSet<>();
	}

	public synchronized void addConnection(int src, int dest, double score) {
		if (!this.adjacents.containsKey(src)) {
			Set<Edge> edges = new TreeSet<>();
			edges.add(new Edge(dest, score));
			this.adjacents.put(src, edges);
		} else {
			Set<Edge> edges = this.adjacents.get(src);
			int sizeBefore = edges.size();
			edges.add(new Edge(dest, score));
			totalEdges += edges.size() - sizeBefore;
			nodes.add(dest);
			totalNodes = nodes.size();
		}
	}

	public synchronized Set<Integer> getNodeList() {
		return nodes;
	}

	public synchronized long countNodes() {
		return totalNodes;
	}

	public synchronized long countEdges() {
		return totalEdges;
	}

	public synchronized PangeneNet getUndirected() {
		PangeneNet pnet = new PangeneNet();

		for (Map.Entry<Integer, Set<Edge>> adjs_e : this.adjacents.entrySet()) {
			for (Edge i : adjs_e.getValue()) {
				pnet.addConnection(adjs_e.getKey(), i.node, i.score);
				pnet.addConnection(i.node, adjs_e.getKey(), i.score);
			}
		}

		return pnet;
	}

	public synchronized List<Set<Integer>> undirectedConnectedComponents() {
		PangeneNet pnet = this.getUndirected();
		return pnet.connectedComponents();
	}

	public synchronized List<Set<Integer>> connectedComponents() {
		List<Set<Integer>> coco = new LinkedList<>();

		Set<Integer> nodes = this.getNodeList();
		HashMap<Integer, Boolean> visited = new HashMap<>();
		for (Integer n : nodes) {
			visited.put(n, false);
		}

		for (Integer root : nodes) {
			if (!visited.get(root)) {
				List<Integer> co = new LinkedList<>();
				co.add(root);
				visited.put(root, true);

				Set<Integer> coconodes = new HashSet<>();
				coconodes.add(root);

				for (int i = 0; i < co.size(); i++) {
					visited.put(co.get(i), true);
					if (this.adjacents.containsKey(co.get(i))) {
						for (Edge ed : this.adjacents.get(co.get(i))) {
							if (!coconodes.contains(ed.node)) {
								co.add(ed.node);
								coconodes.add(ed.node);
							}
							//if(!co.contains(ed.node)){
							//	co.add(ed.node);
							//}
						}
					}
				}

				coco.add(new HashSet<>(co));
			}
		}

		return coco;
	}

	synchronized void printDegreeDistribution() {
		Map<Integer, Integer> distr = new TreeMap<>();
		for (Map.Entry<Integer, Set<Edge>> adjs_e : this.adjacents.entrySet()) {
			distr.put(adjs_e.getValue().size(), distr.getOrDefault(adjs_e.getValue().size(), 0) + 1);
		}
		for (Map.Entry<Integer, Integer> en : distr.entrySet()) {
			System.out.println(en.getKey() + "\t" + en.getValue());
		}
	}

	synchronized void printUndirectedDegreeDistribution() {
		PangeneNet pnet = this.getUndirected();
		pnet.printDegreeDistribution();
	}

	synchronized void print() {
		for (Map.Entry<Integer, Set<Edge>> adjs_e : this.adjacents.entrySet()) {
			System.out.print(adjs_e.getKey() + ": ");
			for (Edge ed : adjs_e.getValue()) {
				System.out.print("(" + ed.node + "," + ed.score + ")");
			}
			System.out.print("\n");
		}
	}

	synchronized void saveToFile(String file, boolean directed) throws Exception {
		BufferedWriter writer = new BufferedWriter(new FileWriter(new File(file)));
		if (directed) {
			for (Map.Entry<Integer, Set<Edge>> adjs_e : this.adjacents.entrySet()) {
				for (Edge ed : adjs_e.getValue()) {
					writer.write(adjs_e.getKey() + "\t" + ed.node + "\t" + ed.score + "\n");
				}
			}
		} else {
			for (Map.Entry<Integer, Set<Edge>> adjs_e : this.adjacents.entrySet()) {
				for (Edge ed : adjs_e.getValue()) {
					if (adjs_e.getKey() <= ed.node) {
						writer.write(adjs_e.getKey() + "\t" + ed.node + "\t" + ed.score + "\n");
					}
				}
			}
		}

		writer.flush();
		writer.close();
	}
}

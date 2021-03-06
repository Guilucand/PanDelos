package infoasys.cli.pangenes;

import java.io.BufferedReader;
import java.io.FileReader;
import java.util.HashMap;
import java.util.Map;
import java.util.Set;
import java.util.TreeSet;
import java.util.Vector;

public class PangeneIData {

    public Vector<String> sequences = new Vector<>();
    public Vector<String> sequenceName = new Vector<>();
    public Vector<String> sequenceDescription = new Vector<>();

    public Vector<Integer> sequenceGenome = new Vector<>();
    public Vector<String> genomeNames = new Vector<>();

    public Set<Character> getAlphabet(){
        Set<Character> alphabet = new TreeSet<>();
        for(String s : this.sequences){
            for(int i=0; i<s.length(); i++){
                alphabet.add(s.charAt(i));
            }
        }
        return alphabet;
    }

    public static PangeneIData readFromFile(String file) throws Exception{
        PangeneIData d = new PangeneIData();

        BufferedReader br = new BufferedReader(new FileReader(file));
        String line;
        boolean nameLine = true;
        String seqName=null, genomeName=null, product=null, seq;
        String[] cc;
        Map<String,Integer> genomeID = new HashMap<>();
        Integer genomeid;
        while ((line = br.readLine()) != null) {

            if (line.trim().isEmpty()) {
                continue;
            }

            // process the line.
            if(nameLine){
                cc = line.trim().split("\t");
                genomeName = cc[0];
                seqName = cc[1];
                product = cc[2];
            }
            else{
                seq = line.trim();
                d.sequences.add(seq);
                d.sequenceName.add(seqName);
                genomeid = genomeID.get(genomeName);
                if(genomeid == null){
                    genomeid = genomeID.size();
                    genomeID.put(genomeName, genomeid);
                }
                d.sequenceGenome.add(genomeid);
                d.sequenceDescription.add(product);
            }
            nameLine = !nameLine;
        }

        d.genomeNames = new Vector<>();
        d.genomeNames.setSize(genomeID.size());
        for(Map.Entry<String,Integer> en : genomeID.entrySet()){
            d.genomeNames.set(en.getValue(), en.getKey());
        }

        return d;
    }

    public Map<String,Integer> getNamedGenomeLengths(){
        Map<String,Integer> lengths = new HashMap<>();
        for(int i=0; i<this.sequences.size(); i++){
            lengths.put( this.genomeNames.get(i),  lengths.getOrDefault(this.genomeNames.get(i), 0) + this.sequences.get(i).length());
        }
        return lengths;
    }
    public Map<Integer,Integer> getGenomeLengths(){
        Map<Integer,Integer> lengths = new HashMap<>();
        for(int i=0; i<this.sequences.size(); i++){
            lengths.put( this.sequenceGenome.get(i),  lengths.getOrDefault(this.sequenceGenome.get(i), 0) + this.sequences.get(i).length());
        }
        return lengths;
    }
    public Map<Integer,Integer> getGenomeKLengths(int k){
        Map<Integer,Integer> lengths = new HashMap<>();
        for(int i=0; i<this.sequences.size(); i++){
            lengths.put( this.sequenceGenome.get(i),  lengths.getOrDefault(this.sequenceGenome.get(i), 0) + this.sequences.get(i).length() -k +1);
        }
        return lengths;
    }

    public Map<Integer, Vector<Integer> > getGenomeSets(){
        Map<Integer, Vector<Integer> > genomes = new HashMap<>();
        for(int g=0; g<genomeNames.size(); g++){
            genomes.put(g, new Vector<>());
        }
        for(int i=0; i<sequenceGenome.size(); i++){
            genomes.get( sequenceGenome.get(i) ).add(i);
        }
        return genomes;
    }
}

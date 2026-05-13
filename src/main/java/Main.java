import java.io.*;
import java.util.*;
import java.lang.Math;

import org.biojava.nbio.core.exceptions.CompoundNotFoundException;
import org.biojava.nbio.core.sequence.AccessionID;
import org.biojava.nbio.core.sequence.io.FastaWriterHelper;
import picocli.CommandLine;
import org.biojava.nbio.core.sequence.DNASequence;
import org.biojava.nbio.core.sequence.io.DNASequenceCreator;
import org.biojava.nbio.core.sequence.compound.AmbiguityDNACompoundSet; 
import org.biojava.nbio.core.sequence.io.FastaReader;
import org.biojava.nbio.core.sequence.io.GenericFastaHeaderParser;

@CommandLine.Command(name = "dunes", mixinStandardHelpOptions = true, version = "0.0")
public class Main implements Runnable{
    @CommandLine.Option(names={"-i", "--inFile"}, description="input fasta file with viral sequences to be mutated",
            paramLabel = "FILE", required=true)
    private File inputFile;
    @CommandLine.Option(names={"-o", "--outFile"},
            description="fasta file with mutated sequences (default: \"<inFile base>\"Mut.fasta\")",
            paramLabel = "<outputFolder>")
    private File outputFile;
    @CommandLine.Option(names={"-m", "--mutation-rate"},
            description="mutation rate in substitutions per nucleotide per year (s/n/y) (default: 4.1e-3)",
            paramLabel = "m")
    private double rate = 4.1e-3;
    @CommandLine.Option(names={"-y", "--years"},
            description="years of evolution (default: 1)",
            paramLabel = "y")
    private double years = 1;
    @CommandLine.Option(names={"-n", "--mutants-number"},
            description="number of mutants for a strain (default: 1)",
            paramLabel = "n")
    private int n = 1;
    @CommandLine.Option(names={"-d", "--distribution"},
            description="mutation model to use: 'simple' (equal-probability) or 'hiv' (HIV-like) (default: 'simple')",
            paramLabel = "distribution")
    private String distribution = "simple";

    private double nuc_mut_prob;
    private static Random rand = new Random();

    private static final List<Character> CANONICAL = Arrays.asList('A', 'C', 'G', 'T');

    // Matrix to use for mutation probabilities
    private double[][] selectedMatrix;

    private static final Map<Character, char[]> IUPAC_MAP = new HashMap<Character, char[]>() {{
        put('A', new char[]{'A'});
        put('C', new char[]{'C'});
        put('G', new char[]{'G'});
        put('T', new char[]{'T'});
        put('U', new char[]{'T'});
        put('R', new char[]{'A','G'});
        put('Y', new char[]{'C','T'});
        put('S', new char[]{'G','C'});
        put('W', new char[]{'A','T'});
        put('K', new char[]{'G','T'});
        put('M', new char[]{'A','C'});
        put('B', new char[]{'C','G','T'});
        put('D', new char[]{'A','G','T'});
        put('H', new char[]{'A','C','T'});
        put('V', new char[]{'A','C','G'});
        put('N', new char[]{'A','C','G','T'});
        put('-', new char[]{'-'});
    }};

    // Equal-probability simple model
    private static final double[][] SIMPLE_MATRIX = {
        {0.0, 1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0}, // A -> C,G,T
        {1.0 / 3.0, 0.0, 1.0 / 3.0, 1.0 / 3.0}, // C -> A,G,T
        {1.0 / 3.0, 1.0 / 3.0, 0.0, 1.0 / 3.0}, // G -> A,C,T
        {1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0, 0.0} // T -> A,C,G
    };

    // First-pass HIV-like matrix
    private static final double[][] HIV_MATRIX = {
        {0.0, 0.15, 0.70, 0.15}, // A -> C,G,T
        {0.15, 0.0, 0.15, 0.70}, // C -> A,G,T
        {0.80, 0.10, 0.0, 0.10}, // G -> A,C,T
        {0.15, 0.70, 0.15, 0.0} // T -> A,C,G
    };


    public void run() {
        try {
            nuc_mut_prob = Math.pow(1+rate,years) - 1;

            // Select mutation matrix based on distribution selected
            if (distribution.equalsIgnoreCase("hiv")) {
                selectedMatrix = HIV_MATRIX;
            } else {
                selectedMatrix = SIMPLE_MATRIX;
            }

            AmbiguityDNACompoundSet ambiguityDNACompoundSet = AmbiguityDNACompoundSet.getDNACompoundSet();
            DNASequenceCreator ambigDNASequenceCreator = new DNASequenceCreator(ambiguityDNACompoundSet);
            GenericFastaHeaderParser fastaHeaderParser = new GenericFastaHeaderParser();
            FastaReader fastaReader = new FastaReader(inputFile, fastaHeaderParser, ambigDNASequenceCreator);
            LinkedHashMap<String, DNASequence> seqs = fastaReader.process(); 
            Set<String> seq_names = seqs.keySet();
            LinkedHashSet<DNASequence> out_seqs = new LinkedHashSet<>();
            if(outputFile == null) {
                String[] tokens = inputFile.toString().split("\\.(?=[^\\.]+$)");
                outputFile = new File(String.format("%sMut.%s",tokens[0],tokens[1]));
            }
            for (String seq_name : seq_names) {
                 out_seqs.addAll(mutate_seqs(seqs.get(seq_name), seq_name));
            }
            FastaWriterHelper.writeNucleotideSequence(outputFile, out_seqs);
        }
        catch(Exception e) {
            e.printStackTrace();
        }
    }

    public static void main(String[] args) {
        CommandLine.run(new Main(), System.out, args);
    }

    private LinkedHashSet<DNASequence> mutate_seqs(DNASequence seq, String seq_name)
            throws CompoundNotFoundException {
        LinkedHashSet<DNASequence> out_seqs = new LinkedHashSet<>();
        for(int i=0; i<n; i++) {
            DNASequence mut_seq = mutate_seq(seq);
            mut_seq.setAccession(new AccessionID(String.format("%s_%d",seq_name,i)));
            out_seqs.add(mut_seq);
        }
        return out_seqs;
    }

    private static char resolveAmbiguousBase(char nuc) {
        nuc = Character.toUpperCase(nuc);
        char[] choices = IUPAC_MAP.get(nuc);

        if (choices == null) {
            throw new IllegalArgumentException("Unsupported nucleotide code: " + nuc);
        }

        return choices[rand.nextInt(choices.length)];
    }

    private char mutateCanonicalBase(char base) {
        base = Character.toUpperCase(base);
        int baseIdx = CANONICAL.indexOf(base);
        if (baseIdx == -1) {
            throw new IllegalArgumentException("Base must be one of A, C, G, T");
        }
        double[] probs = selectedMatrix[baseIdx];
        // Build cumulative distribution
        double r = rand.nextDouble();
        double cumulative = 0.0;
        for (int i = 0; i < probs.length; i++) {
            if (i == baseIdx) continue; // skip self
            cumulative += probs[i];
            if (r < cumulative) {
                return CANONICAL.get(i);
            }
        }
        // Fallback (should not happen)
        for (int i = 0; i < probs.length; i++) {
            if (i != baseIdx && probs[i] > 0) {
                return CANONICAL.get(i);
            }
        }
        throw new IllegalStateException("No valid mutation found for base: " + base);
    }

    private DNASequence mutate_seq(DNASequence seq) throws CompoundNotFoundException {
        StringBuilder mut_seq = new StringBuilder();
        for (int i=0; i<seq.getLength(); i++) {
            char nuc = Character.toUpperCase(seq.getCompoundAt(i+1).toString().charAt(0));

            // Preserve gaps unchanged
            if (nuc == '-') {
                mut_seq.append(nuc);
                continue;
            }

            // No mutation event: keep original symbol exactly as written
            if (Math.random() > nuc_mut_prob) {
                mut_seq.append(nuc);
                continue;
            }

            // Mutation event: resolve ambiguity to one compatible base, then mutate away from that base
            char resolved = resolveAmbiguousBase(nuc);
            char mutated = mutateCanonicalBase(resolved);
            mut_seq.append(mutated);
        }
        return new DNASequence(mut_seq.toString(), AmbiguityDNACompoundSet.getDNACompoundSet());
    }
}

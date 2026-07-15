import org.biojava.nbio.core.exceptions.CompoundNotFoundException;
import org.biojava.nbio.core.sequence.AccessionID;
import org.biojava.nbio.core.sequence.DNASequence;
import org.biojava.nbio.core.sequence.compound.AmbiguityDNACompoundSet;
import org.biojava.nbio.core.sequence.io.DNASequenceCreator;
import org.biojava.nbio.core.sequence.io.FastaReader;
import org.biojava.nbio.core.sequence.io.FastaWriterHelper;
import org.biojava.nbio.core.sequence.io.GenericFastaHeaderParser;

import java.io.File;
import java.io.FileInputStream;
import java.io.InputStream;
import java.io.PrintWriter;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.HashMap;
import java.util.Random;

public class DunesEngine {

    public interface ProgressListener {
        void onProgress(String message, int currentParent, int totalParents);
    }

    public static class QCResult {
        public final String parentSequence;
        public final String descendantSequence;
        public final int sequenceLength;
        public final int rawDifferences;
        public final double normalizedDistance;
        public final String metric;
        public final double threshold;
        public final boolean pass;

        public QCResult(String parentSequence,
                        String descendantSequence,
                        int sequenceLength,
                        int rawDifferences,
                        double normalizedDistance,
                        String metric,
                        double threshold,
                        boolean pass) {
            this.parentSequence = parentSequence;
            this.descendantSequence = descendantSequence;
            this.sequenceLength = sequenceLength;
            this.rawDifferences = rawDifferences;
            this.normalizedDistance = normalizedDistance;
            this.metric = metric;
            this.threshold = threshold;
            this.pass = pass;
        }
    }

    private static final Random rand = new Random();

    private static final Map<Character, Integer> BASE_INDEX = new HashMap<Character, Integer>() {{
        put('A', 0);
        put('C', 1);
        put('G', 2);
        put('T', 3);
    }};

    private static final Map<Character, char[]> IUPAC_MAP = new HashMap<Character, char[]>() {{
        put('A', new char[]{'A'});
        put('C', new char[]{'C'});
        put('G', new char[]{'G'});
        put('T', new char[]{'T'});
        put('U', new char[]{'T'});
        put('R', new char[]{'A', 'G'});
        put('Y', new char[]{'C', 'T'});
        put('S', new char[]{'G', 'C'});
        put('W', new char[]{'A', 'T'});
        put('K', new char[]{'G', 'T'});
        put('M', new char[]{'A', 'C'});
        put('B', new char[]{'C', 'G', 'T'});
        put('D', new char[]{'A', 'G', 'T'});
        put('H', new char[]{'A', 'C', 'T'});
        put('V', new char[]{'A', 'C', 'G'});
        put('N', new char[]{'A', 'C', 'G', 'T'});
        put('-', new char[]{'-'});
    }};

    private static final double[][] SIMPLE_MATRIX = {
            {0.0, 1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0},
            {1.0 / 3.0, 0.0, 1.0 / 3.0, 1.0 / 3.0},
            {1.0 / 3.0, 1.0 / 3.0, 0.0, 1.0 / 3.0},
            {1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0, 0.0}
    };

    private static final double[][] HIV_MATRIX = {
            {0.0, 0.15, 0.70, 0.15},
            {0.15, 0.0, 0.15, 0.70},
            {0.80, 0.10, 0.0, 0.10},
            {0.15, 0.70, 0.15, 0.0}
    };

    private static final double[][] HIV_APOBEC_MATRIX = {
            {0.0, 0.15, 0.70, 0.15},
            {0.15, 0.0, 0.15, 0.70},
            {0.95, 0.03, 0.0, 0.02},
            {0.15, 0.70, 0.15, 0.0}
    };

    private static final Map<String, Character> CODON_TABLE = new HashMap<String, Character>() {{
        put("TTT", 'F'); put("TTC", 'F'); put("TTA", 'L'); put("TTG", 'L');
        put("TCT", 'S'); put("TCC", 'S'); put("TCA", 'S'); put("TCG", 'S');
        put("TAT", 'Y'); put("TAC", 'Y'); put("TAA", '*'); put("TAG", '*');
        put("TGT", 'C'); put("TGC", 'C'); put("TGA", '*'); put("TGG", 'W');

        put("CTT", 'L'); put("CTC", 'L'); put("CTA", 'L'); put("CTG", 'L');
        put("CCT", 'P'); put("CCC", 'P'); put("CCA", 'P'); put("CCG", 'P');
        put("CAT", 'H'); put("CAC", 'H'); put("CAA", 'Q'); put("CAG", 'Q');
        put("CGT", 'R'); put("CGC", 'R'); put("CGA", 'R'); put("CGG", 'R');

        put("ATT", 'I'); put("ATC", 'I'); put("ATA", 'I'); put("ATG", 'M');
        put("ACT", 'T'); put("ACC", 'T'); put("ACA", 'T'); put("ACG", 'T');
        put("AAT", 'N'); put("AAC", 'N'); put("AAA", 'K'); put("AAG", 'K');
        put("AGT", 'S'); put("AGC", 'S'); put("AGA", 'R'); put("AGG", 'R');

        put("GTT", 'V'); put("GTC", 'V'); put("GTA", 'V'); put("GTG", 'V');
        put("GCT", 'A'); put("GCC", 'A'); put("GCA", 'A'); put("GCG", 'A');
        put("GAT", 'D'); put("GAC", 'D'); put("GAA", 'E'); put("GAG", 'E');
        put("GGT", 'G'); put("GGC", 'G'); put("GGA", 'G'); put("GGG", 'G');
    }};

    private int lastParentCount = 0;
    private int lastDescendantCount = 0;
    private final List<QCResult> lastQcResults = new ArrayList<QCResult>();
    private int lastQcPassCount = 0;
    private int lastQcFailCount = 0;

    public void run(File inFile,
                    File outFile,
                    double mutationRate,
                    double years,
                    int n,
                    String model,
                    boolean codonAware,
                    boolean apobec,
                    double apobecRate,
                    File siteWeightsFile) throws Exception {
        run(inFile, outFile, mutationRate, years, n, model, codonAware, apobec, apobecRate,
                siteWeightsFile, null, false, "normalized_nt_distance", 0.0);
    }

    public void run(File inFile,
                    File outFile,
                    double mutationRate,
                    double years,
                    int n,
                    String model,
                    boolean codonAware,
                    boolean apobec,
                    double apobecRate,
                    File siteWeightsFile,
                    ProgressListener progressListener,
                    boolean qcEnabled,
                    String qcMetric,
                    double qcThreshold) throws Exception {

        if (n < 1) {
            throw new IllegalArgumentException("Mutants per sequence (-n) must be >= 1.");
        }
        if (mutationRate < 0.0) {
            throw new IllegalArgumentException("Mutation rate must be >= 0.");
        }
        if (years < 0.0) {
            throw new IllegalArgumentException("Years must be >= 0.");
        }
        if (apobecRate < 0.0 || apobecRate > 1.0) {
            throw new IllegalArgumentException("APOBEC rate must be between 0 and 1.");
        }
        if (!"simple".equalsIgnoreCase(model) && !"hiv".equalsIgnoreCase(model)) {
            throw new IllegalArgumentException("Model must be 'simple' or 'hiv'.");
        }
        if (!"normalized_nt_distance".equals(qcMetric) && !"raw_nt_differences".equals(qcMetric)) {
            throw new IllegalArgumentException("QC metric must be 'normalized_nt_distance' or 'raw_nt_differences'.");
        }

        double nucMutProb = Math.pow(1.0 + mutationRate, years) - 1.0;
        if (nucMutProb > 1.0) {
            nucMutProb = 1.0;
        }

        Map<String, DNASequence> seqs;
        try (InputStream is = new FileInputStream(inFile)) {
            FastaReader<DNASequence, org.biojava.nbio.core.sequence.compound.NucleotideCompound> fastaReader =
                    new FastaReader<DNASequence, org.biojava.nbio.core.sequence.compound.NucleotideCompound>(
                            is,
                            new GenericFastaHeaderParser<DNASequence, org.biojava.nbio.core.sequence.compound.NucleotideCompound>(),
                            new DNASequenceCreator(AmbiguityDNACompoundSet.getDNACompoundSet())
                    );
            seqs = fastaReader.process();
        }

        lastParentCount = seqs.size();
        lastDescendantCount = 0;
        lastQcResults.clear();
        lastQcPassCount = 0;
        lastQcFailCount = 0;

        List<DNASequence> outSeqs = new ArrayList<DNASequence>();
        int totalParents = seqs.size();
        int parentIndex = 0;

        for (Map.Entry<String, DNASequence> entry : seqs.entrySet()) {
            if (Thread.currentThread().isInterrupted()) {
                throw new InterruptedException("DUNES run cancelled.");
            }

            parentIndex++;
            String seqName = entry.getKey();

            if (progressListener != null) {
                progressListener.onProgress(
                        "Processing parent " + parentIndex + " of " + totalParents + ": " + seqName,
                        parentIndex,
                        totalParents
                );
            }

            DNASequence seq = entry.getValue();
            double[] siteWeights = loadSiteWeights(siteWeightsFile, seq.getLength());
            List<DNASequence> generated = mutateSeqs(
                    seq, seqName, n, nucMutProb, model, codonAware, apobec, apobecRate,
                    siteWeights, qcEnabled, qcMetric, qcThreshold
            );
            outSeqs.addAll(generated);
            lastDescendantCount += generated.size();
        }

        FastaWriterHelper.writeNucleotideSequence(outFile, outSeqs);
    }

    private List<DNASequence> mutateSeqs(DNASequence seq,
                                         String seqName,
                                         int n,
                                         double nucMutProb,
                                         String model,
                                         boolean codonAware,
                                         boolean apobec,
                                         double apobecRate,
                                         double[] siteWeights,
                                         boolean qcEnabled,
                                         String qcMetric,
                                         double qcThreshold) throws CompoundNotFoundException, InterruptedException {

        List<DNASequence> outSeqs = new ArrayList<DNASequence>();
        String parentSeqString = seq.getSequenceAsString().toUpperCase();

        for (int i = 0; i < n; i++) {
            if (Thread.currentThread().isInterrupted()) {
                throw new InterruptedException("DUNES run cancelled.");
            }

            boolean apobecMode = apobec && rand.nextDouble() < apobecRate;
            DNASequence mutSeq = mutateSeq(seq, nucMutProb, model, codonAware, apobecMode, siteWeights);
            mutSeq.setAccession(new AccessionID(String.format("%s_%d", seqName, i)));
            outSeqs.add(mutSeq);

            if (qcEnabled) {
                String childSeqString = mutSeq.getSequenceAsString().toUpperCase();
                int rawDiffs = countDifferences(parentSeqString, childSeqString);
                double normDist = normalizedDistance(parentSeqString, childSeqString);

                boolean pass;
                if ("raw_nt_differences".equals(qcMetric)) {
                    pass = rawDiffs <= qcThreshold;
                } else {
                    pass = normDist <= qcThreshold;
                }

                if (pass) {
                    lastQcPassCount++;
                } else {
                    lastQcFailCount++;
                }

                lastQcResults.add(new QCResult(
                        seqName,
                        mutSeq.getAccession().toString(),
                        childSeqString.length(),
                        rawDiffs,
                        normDist,
                        qcMetric,
                        qcThreshold,
                        pass
                ));
            }
        }

        return outSeqs;
    }

    private double[] loadSiteWeights(File f, int expectedLen) throws Exception {
        double[] w = new double[expectedLen];
        Arrays.fill(w, 1.0);

        if (f == null) {
            return w;
        }

        List<String> lines = Files.readAllLines(f.toPath());
        if (lines.size() != expectedLen) {
            throw new IllegalArgumentException("Site weights length (" + lines.size()
                    + ") does not match sequence length (" + expectedLen + ").");
        }

        for (int i = 0; i < expectedLen; i++) {
            String s = lines.get(i).trim();
            if (s.isEmpty()) {
                throw new IllegalArgumentException("Blank line in site weights file at line " + (i + 1));
            }
            w[i] = Double.parseDouble(s);
            if (w[i] < 0.0) {
                w[i] = 0.0;
            }
        }

        return w;
    }

    private DNASequence mutateSeq(DNASequence seq,
                                  double nucMutProb,
                                  String model,
                                  boolean codonAware,
                                  boolean apobecMode,
                                  double[] siteWeights) throws CompoundNotFoundException, InterruptedException {

        char[] working = seq.getSequenceAsString().toUpperCase().toCharArray();

        for (int i = 0; i < working.length; i++) {
            if (Thread.currentThread().isInterrupted()) {
                throw new InterruptedException("DUNES run cancelled.");
            }

            char nuc = working[i];

            if (nuc == '-') {
                continue;
            }

            double effectiveProb = nucMutProb * siteWeights[i];
            if (effectiveProb > 1.0) {
                effectiveProb = 1.0;
            }

            if (rand.nextDouble() > effectiveProb) {
                continue;
            }

            char resolved = resolveAmbiguousBase(nuc);
            char mutated = chooseMutationUsingMatrix(resolved, model, apobecMode);

            if (!codonAware || acceptMutationByCodonContext(working, i, mutated)) {
                working[i] = mutated;
            }
        }

        return new DNASequence(new String(working), AmbiguityDNACompoundSet.getDNACompoundSet());
    }

    private static char resolveAmbiguousBase(char nuc) {
        nuc = Character.toUpperCase(nuc);
        char[] choices = IUPAC_MAP.get(nuc);

        if (choices == null) {
            throw new IllegalArgumentException("Unsupported nucleotide code: " + nuc);
        }

        return choices[rand.nextInt(choices.length)];
    }

    private char chooseMutationUsingMatrix(char base, String model, boolean apobecMode) {
        base = Character.toUpperCase(base);

        double[][] matrix = SIMPLE_MATRIX;
        if ("hiv".equalsIgnoreCase(model)) {
            matrix = apobecMode ? HIV_APOBEC_MATRIX : HIV_MATRIX;
        }

        Integer idx = BASE_INDEX.get(base);
        if (idx == null) {
            throw new IllegalArgumentException("Unsupported canonical base: " + base);
        }

        double r = rand.nextDouble();
        double cum = 0.0;
        char[] bases = {'A', 'C', 'G', 'T'};

        for (int j = 0; j < 4; j++) {
            cum += matrix[idx][j];
            if (r <= cum) {
                return bases[j];
            }
        }

        for (int j = 0; j < 4; j++) {
            if (j != idx) {
                return bases[j];
            }
        }

        return base;
    }

    private boolean acceptMutationByCodonContext(char[] seqChars, int pos, char newBase) {
        int codonStart = (pos / 3) * 3;
        if (codonStart + 2 >= seqChars.length) {
            return true;
        }

        char[] oldCodon = new char[]{
                Character.toUpperCase(seqChars[codonStart]),
                Character.toUpperCase(seqChars[codonStart + 1]),
                Character.toUpperCase(seqChars[codonStart + 2])
        };

        if (!isCanonicalCodon(oldCodon)) {
            return true;
        }

        char[] newCodon = oldCodon.clone();
        newCodon[pos - codonStart] = Character.toUpperCase(newBase);

        char aaOld = translateCodon(oldCodon);
        char aaNew = translateCodon(newCodon);

        if (aaOld == '?' || aaNew == '?') {
            return true;
        }
        if (aaNew == '*') {
            return false;
        }
        if (aaOld == aaNew) {
            return true;
        }

        double p = isConservativeChange(aaOld, aaNew) ? 0.5 : 0.2;
        return rand.nextDouble() < p;
    }

    private boolean isCanonicalCodon(char[] codon) {
        if (codon.length != 3) {
            return false;
        }
        for (char c : codon) {
            if (!BASE_INDEX.containsKey(Character.toUpperCase(c))) {
                return false;
            }
        }
        return true;
    }

    private char translateCodon(char[] codon) {
        String key = new String(codon).toUpperCase();
        Character aa = CODON_TABLE.get(key);
        return aa == null ? '?' : aa;
    }

    private boolean isConservativeChange(char aaOld, char aaNew) {
        return aminoAcidGroup(aaOld) == aminoAcidGroup(aaNew);
    }

    private int aminoAcidGroup(char aa) {
        aa = Character.toUpperCase(aa);

        if ("AVLIMFWY".indexOf(aa) >= 0) return 1;
        if ("STNQCGP".indexOf(aa) >= 0) return 2;
        if ("KRH".indexOf(aa) >= 0) return 3;
        if ("DE".indexOf(aa) >= 0) return 4;
        return 0;
    }

    private int countDifferences(String parent, String child) {
        int len = Math.min(parent.length(), child.length());
        int diffs = Math.abs(parent.length() - child.length());

        for (int i = 0; i < len; i++) {
            if (Character.toUpperCase(parent.charAt(i)) != Character.toUpperCase(child.charAt(i))) {
                diffs++;
            }
        }
        return diffs;
    }

    private double normalizedDistance(String parent, String child) {
        int maxLen = Math.max(parent.length(), child.length());
        if (maxLen == 0) {
            return 0.0;
        }
        return countDifferences(parent, child) / (double) maxLen;
    }

    public void writeLastQcResultsCsv(File csvFile) throws Exception {
        try (PrintWriter pw = new PrintWriter(csvFile)) {
            pw.println("parent_sequence,descendant_sequence,sequence_length,raw_differences,normalized_distance,metric,threshold,pass");

            for (QCResult r : lastQcResults) {
                pw.println(
                        escapeCsv(r.parentSequence) + "," +
                        escapeCsv(r.descendantSequence) + "," +
                        r.sequenceLength + "," +
                        r.rawDifferences + "," +
                        r.normalizedDistance + "," +
                        escapeCsv(r.metric) + "," +
                        r.threshold + "," +
                        (r.pass ? "Yes" : "No")
                );
            }
        }
    }

    private String escapeCsv(String s) {
        if (s == null) {
            return "";
        }
        if (s.contains(",") || s.contains("\"") || s.contains("\n")) {
            return "\"" + s.replace("\"", "\"\"") + "\"";
        }
        return s;
    }

    public static File defaultQcCsvPath(File outFile) {
        String outPath = outFile.getAbsolutePath();
        if (outPath.toLowerCase().endsWith(".fasta")) {
            return new File(outPath.substring(0, outPath.length() - 6) + "_QC.csv");
        }
        if (outPath.toLowerCase().endsWith(".fa")) {
            return new File(outPath.substring(0, outPath.length() - 3) + "_QC.csv");
        }
        return new File(outPath + "_QC.csv");
    }

    public int getLastParentCount() {
        return lastParentCount;
    }

    public int getLastDescendantCount() {
        return lastDescendantCount;
    }

    public List<QCResult> getLastQcResults() {
        return new ArrayList<QCResult>(lastQcResults);
    }

    public int getLastQcPassCount() {
        return lastQcPassCount;
    }

    public int getLastQcFailCount() {
        return lastQcFailCount;
    }
}

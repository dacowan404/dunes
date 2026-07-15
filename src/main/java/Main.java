import picocli.CommandLine.Command;
import picocli.CommandLine.Option;
import picocli.CommandLine;

import java.io.File;

@Command(
        name = "dunes",
        mixinStandardHelpOptions = true,
        version = "0.9",
        description = "DUNES: sequence mutation simulator"
)
public class Main implements Runnable {

    public static void main(String[] args) {
        CommandLine.run(new Main(), args);
    }

    @Option(names = {"-i", "--inFile"}, required = true, description = "Input FASTA file")
    private File inFile;

    @Option(names = {"-o", "--outFile"}, description = "Output FASTA file")
    private File outFile;

    @Option(names = {"-m", "--mutation-rate"}, description = "Mutation rate per nucleotide per year")
    private double mutationRate = 0.004;

    @Option(names = {"-y", "--years"}, description = "Years of evolution")
    private double years = 2.0;

    @Option(names = {"-n", "--mutants-number"}, description = "Mutants per input sequence")
    private int n = 10;

    @Option(names = {"--model"}, description = "Evolution model: simple or HIV")
    private String model = "hiv";

    @Option(names = {"--codon-aware"}, description = "Use codon-aware acceptance filter")
    private boolean codonAware = false;

    @Option(names = {"--apobec"}, description = "Enable APOBEC-like G->A hypermutation")
    private boolean apobec = false;

    @Option(names = {"--apobec-rate"}, description = "Probability a descendant is APOBEC-affected")
    private double apobecRate = 0.02;

    @Option(names = {"--site-weights"}, description = "Optional per-site weights file")
    private File siteWeightsFile;

    @Option(names = {"--qc-enable"}, description = "Enable QC distance check")
    private boolean qcEnabled = false;

    @Option(names = {"--qc-metric"}, description = "QC metric: normalized_nt_distance or raw_nt_differences")
    private String qcMetric = "normalized_nt_distance";

    @Option(names = {"--qc-threshold"}, description = "Maximum allowed QC threshold")
    private double qcThreshold = 0.01;

    @Option(names = {"--qc-write-csv"}, description = "Write QC CSV summary next to output FASTA")
    private boolean qcWriteCsv = true;

    @Override
    public void run() {
        try {
            if (outFile == null) {
                String inName = inFile.getName();
                int dot = inName.lastIndexOf('.');
                String base = (dot > 0) ? inName.substring(0, dot) : inName;

                File parent = inFile.getAbsoluteFile().getParentFile();
                if (parent == null) {
                    parent = new File(".");
                }

                outFile = new File(parent, base + "Mut.fasta");
            }

            String normalizedModel = "HIV".equalsIgnoreCase(model) ? "hiv" : model.toLowerCase();

            DunesEngine engine = new DunesEngine();
            engine.run(
                    inFile,
                    outFile,
                    mutationRate,
                    years,
                    n,
                    normalizedModel,
                    codonAware,
                    apobec,
                    apobecRate,
                    siteWeightsFile,
                    null,
                    qcEnabled,
                    qcMetric,
                    qcThreshold
            );

            if (qcEnabled && qcWriteCsv) {
                engine.writeLastQcResultsCsv(DunesEngine.defaultQcCsvPath(outFile));
            }

            System.out.println("Done. Output written to: " + outFile.getAbsolutePath());
            System.out.println("Summary: processed " + engine.getLastParentCount()
                    + " parent sequence(s) and wrote "
                    + engine.getLastDescendantCount() + " descendant sequence(s).");

            if (qcEnabled) {
                System.out.println("QC summary: " + engine.getLastQcPassCount()
                        + " pass, " + engine.getLastQcFailCount() + " fail.");
            }
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
    }
}

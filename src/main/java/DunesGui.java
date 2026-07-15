import javax.swing.BorderFactory;
import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JComboBox;
import javax.swing.JComponent;
import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.ImageIcon;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JProgressBar;
import javax.swing.JScrollPane;
import javax.swing.JTextArea;
import javax.swing.JTextField;
import javax.swing.SwingUtilities;
import javax.swing.SwingWorker;
import javax.swing.ToolTipManager;
import java.awt.BorderLayout;
import java.awt.Desktop;
import java.awt.Image;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Insets;
import java.awt.event.FocusAdapter;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.awt.event.FocusEvent;
import java.io.File;
import java.io.FileWriter;
import java.time.LocalDateTime;
import java.time.format.DateTimeFormatter;

public class DunesGui extends JFrame {

    private static final DateTimeFormatter TS_FORMAT =
            DateTimeFormatter.ofPattern("yyyy-MM-dd HH:mm:ss");
    private static final long MAX_SUGGESTED_FILE_SIZE = 100L * 1024L * 1024L; // 100 MB
    private static final int MIN_RECOMMENDED_JAVA_MAJOR = 11;

    private final JTextField inputField = new JTextField(34);
    private final JTextField outputField = new JTextField(34);
    private final JTextField mutationField = new JTextField("0.004", 10);
    private final JTextField yearsField = new JTextField("2", 10);
    private final JTextField mutantsField = new JTextField("10", 10);
    private final JComboBox<String> modelBox = new JComboBox<String>(new String[]{"simple", "HIV"});
    private final JCheckBox codonAwareBox = new JCheckBox("Codon-aware");
    private final JCheckBox apobecBox = new JCheckBox("APOBEC");
    private final JTextField apobecRateField = new JTextField("0.02", 10);
    private final JTextField siteWeightsField = new JTextField(34);

    private final JCheckBox qcCheckBox = new JCheckBox("Enable QC distance check");
    private final JComboBox<String> qcMetricBox = new JComboBox<String>(new String[]{"normalized_nt_distance", "raw_nt_differences"});
    private final JTextField qcThresholdField = new JTextField("0.01", 10);
    private final JCheckBox qcWriteCsvBox = new JCheckBox("Write QC CSV", true);

    private File lastDirectory = null;

    private final JTextArea logArea = new JTextArea(14, 70);
    private final JButton runButton = new JButton("Run");
    private final JButton cancelButton = new JButton("Cancel");
    private final JButton saveLogButton = new JButton("Save Log");
    private final JButton clearLogButton = new JButton("Clear Log");
    private final JButton resetFormButton = new JButton("Reset Form");
    private final JButton helpButton = new JButton("Help / About");
    private final JButton openManualButton = new JButton("Open Manual");
    private final JProgressBar progressBar = new JProgressBar();

    private SwingWorker<Void, Void> currentWorker = null;

    public DunesGui() {
        super("DUNES");
        setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        setLayout(new BorderLayout(8, 8));

        JPanel form = new JPanel(new GridBagLayout());
        form.setBorder(BorderFactory.createEmptyBorder(10, 10, 10, 10));

        GridBagConstraints c = new GridBagConstraints();
        c.insets = new Insets(5, 5, 5, 5);
        c.fill = GridBagConstraints.HORIZONTAL;

        JButton inputBrowse = new JButton("Browse...");
        inputBrowse.addActionListener(e -> chooseFile(inputField, false, false));

        JButton outputBrowse = new JButton("Browse...");
        outputBrowse.addActionListener(e -> chooseFile(outputField, true, false));

        JButton weightsBrowse = new JButton("Browse...");
        weightsBrowse.addActionListener(e -> chooseFile(siteWeightsField, false, true));

        int row = 0;
        addRow(form, c, row++, "Input FASTA:", inputField, inputBrowse);
        addRow(form, c, row++, "Output FASTA (optional):", outputField, outputBrowse);
        addRow(form, c, row++, "Mutation rate:", mutationField, null);
        addRow(form, c, row++, "Years:", yearsField, null);
        addRow(form, c, row++, "Mutants / sequence:", mutantsField, null);
        addRow(form, c, row++, "Model:", modelBox, null);
        addRow(form, c, row++, "Site weights (optional):", siteWeightsField, weightsBrowse);
        addRow(form, c, row++, "QC metric (optional):", qcMetricBox, null);
        addRow(form, c, row++, "QC threshold (optional):", qcThresholdField, null);

        c.gridx = 0;
        c.gridy = row;
        c.gridwidth = 1;
        c.weightx = 0;
        form.add(new JLabel("Options:"), c);

        JPanel optionsPanel = new JPanel();
        optionsPanel.add(codonAwareBox);
        optionsPanel.add(apobecBox);

        c.gridx = 1;
        c.gridy = row;
        c.gridwidth = 2;
        c.weightx = 1.0;
        form.add(optionsPanel, c);
        row++;

        c.gridx = 0;
        c.gridy = row;
        c.gridwidth = 1;
        c.weightx = 0;
        form.add(new JLabel("QC options:"), c);

        JPanel qcPanel = new JPanel();
        qcPanel.add(qcCheckBox);
        qcPanel.add(qcWriteCsvBox);

        c.gridx = 1;
        c.gridy = row;
        c.gridwidth = 2;
        c.weightx = 1.0;
        form.add(qcPanel, c);
        row++;

        addRow(form, c, row++, "APOBEC rate (optional):", apobecRateField, null);

        runButton.addActionListener(e -> runDunes());

        cancelButton.setEnabled(false);
        cancelButton.addActionListener(e -> cancelRun());

        saveLogButton.addActionListener(e -> saveLog());
        clearLogButton.addActionListener(e -> clearLog());
        resetFormButton.addActionListener(e -> resetForm());
        helpButton.addActionListener(e -> showHelpDialog());
        openManualButton.addActionListener(e -> openManual());

        JPanel buttonPanel = new JPanel();
        buttonPanel.add(runButton);
        buttonPanel.add(cancelButton);
        buttonPanel.add(saveLogButton);
        buttonPanel.add(clearLogButton);
        buttonPanel.add(resetFormButton);
        buttonPanel.add(helpButton);
        buttonPanel.add(openManualButton);

        progressBar.setMinimum(0);
        progressBar.setMaximum(100);
        progressBar.setValue(0);
        progressBar.setStringPainted(true);
        progressBar.setString("0%");

        JPanel southPanel = new JPanel(new BorderLayout());
        southPanel.add(buttonPanel, BorderLayout.NORTH);
        southPanel.add(progressBar, BorderLayout.SOUTH);

        JPanel top = new JPanel(new BorderLayout());

        JLabel logoLabel = createLogoLabel();
        if (logoLabel != null) {
            top.add(logoLabel, BorderLayout.NORTH);
        }

        top.add(form, BorderLayout.CENTER);
        top.add(southPanel, BorderLayout.SOUTH);

        Image windowIcon = loadLogoImage();
        if (windowIcon != null) {
            setIconImage(windowIcon);
        }

        logArea.setEditable(false);
        JScrollPane scrollPane = new JScrollPane(logArea);
        scrollPane.setBorder(BorderFactory.createTitledBorder("Run log"));

        add(top, BorderLayout.NORTH);
        add(scrollPane, BorderLayout.CENTER);

        pack();
        setLocationRelativeTo(null);

        outputField.setToolTipText("Optional. Leave blank to use <input base>Mut.fasta next to the input FASTA.");
        siteWeightsField.setToolTipText("Optional. Leave blank to use equal mutation weight at every site.");
        apobecRateField.setToolTipText("Optional unless APOBEC is enabled. Default: 0.02.");
        qcMetricBox.setToolTipText("Optional unless QC is enabled. Default: normalized_nt_distance.");
        qcThresholdField.setToolTipText("Optional unless QC is enabled. Default: 0.01.");
        configureTooltipBehavior();
        applyFormTooltips(form);

        modelBox.setSelectedItem("HIV");
        qcMetricBox.setSelectedItem("normalized_nt_distance");
        apobecBox.addActionListener(e -> updateApobecFieldState());
        qcCheckBox.addActionListener(e -> updateQcFieldState());
        inputField.addFocusListener(new FocusAdapter() {
            @Override
            public void focusLost(FocusEvent e) {
                updateDefaultOutputFromInput();
            }
        });
        updateApobecFieldState();
        updateQcFieldState();
        showJavaVersionWarningIfNeeded();
    }

    private void showJavaVersionWarningIfNeeded() {
        String javaVersion = System.getProperty("java.version", "unknown");
        int javaMajor = parseJavaMajorVersion(javaVersion);

        if (javaMajor >= 0 && javaMajor < MIN_RECOMMENDED_JAVA_MAJOR) {
            final String message =
                    "DUNES is running with Java " + javaVersion + ".\n\n" +
                    "Java 11 or newer is recommended. Some GUI features may not work correctly on older Java versions.";

            logArea.append("WARNING: DUNES is running with Java " + javaVersion
                    + "; Java 11 or newer is recommended.\n");

            SwingUtilities.invokeLater(() -> JOptionPane.showMessageDialog(
                    DunesGui.this,
                    message,
                    "Java version warning",
                    JOptionPane.WARNING_MESSAGE
            ));
        }
    }

    private int parseJavaMajorVersion(String version) {
        if (version == null) {
            return -1;
        }

        String trimmed = version.trim();
        if (trimmed.startsWith("1.")) {
            return parseLeadingInt(trimmed.substring(2));
        }
        return parseLeadingInt(trimmed);
    }

    private int parseLeadingInt(String value) {
        StringBuilder digits = new StringBuilder();

        for (int i = 0; i < value.length(); i++) {
            char ch = value.charAt(i);
            if (!Character.isDigit(ch)) {
                break;
            }
            digits.append(ch);
        }

        if (digits.length() == 0) {
            return -1;
        }

        try {
            return Integer.parseInt(digits.toString());
        } catch (NumberFormatException ex) {
            return -1;
        }
    }


    private File getApplicationDirectory() {
        try {
            File appLocation = new File(
                    DunesGui.class.getProtectionDomain().getCodeSource().getLocation().toURI()
            ).getCanonicalFile();

            File appDir = appLocation.isFile() ? appLocation.getParentFile() : appLocation;
            if (appDir != null) {
                return appDir;
            }
        } catch (Exception ignored) {
        }

        try {
            return new File(".").getCanonicalFile();
        } catch (Exception ex) {
            return new File(".");
        }
    }

    private Image loadLogoImage() {
        File appDir = getApplicationDirectory();

        File[] candidates = new File[] {
                new File(appDir, "DUNES_logo.png"),
                new File(appDir, "DUNES logo.png"),
                new File(appDir, "DUNES_Logo.png"),
                new File(appDir, "DUNES Logo.png"),
                new File(appDir, "dunes_logo.png")
        };

        for (File f : candidates) {
            if (f.exists()) {
                return new ImageIcon(f.getAbsolutePath()).getImage();
            }
        }

        java.net.URL resource = DunesGui.class.getResource("/DUNES_logo.png");
        if (resource != null) {
            return new ImageIcon(resource).getImage();
        }

        return null;
    }

    private JLabel createLogoLabel() {
        Image logoImage = loadLogoImage();
        if (logoImage == null) {
            return null;
        }

        Image scaled = logoImage.getScaledInstance(180, -1, Image.SCALE_SMOOTH);
        JLabel label = new JLabel(new ImageIcon(scaled));
        label.setBorder(BorderFactory.createEmptyBorder(0, 0, 10, 0));
        label.setHorizontalAlignment(JLabel.CENTER);
        return label;
    }


    private void forceTooltipSupport(JComponent component) {
        if (component == null) {
            return;
        }
        ToolTipManager.sharedInstance().registerComponent(component);
        component.addMouseListener(new MouseAdapter() {
            @Override
            public void mouseEntered(MouseEvent e) {
                ToolTipManager.sharedInstance().mouseMoved(
                        new MouseEvent(component,
                                MouseEvent.MOUSE_MOVED,
                                System.currentTimeMillis(),
                                0,
                                Math.max(1, e.getX()),
                                Math.max(1, e.getY()),
                                0,
                                false)
                );
            }
        });
    }

    private void configureTooltipBehavior() {
        ToolTipManager.sharedInstance().setInitialDelay(250);
        ToolTipManager.sharedInstance().setDismissDelay(15000);
        ToolTipManager.sharedInstance().setReshowDelay(100);
    }

    private void applyFormTooltips(JPanel form) {
        codonAwareBox.setToolTipText("Optional. Uses a simple coding-sequence acceptance filter to reduce disruptive changes.");
        apobecBox.setToolTipText("Optional. Enables APOBEC-like G→A hypermutation.");
        qcCheckBox.setToolTipText("Optional. Enables parent-to-descendant QC distance checking.");
        qcWriteCsvBox.setToolTipText("Optional when QC is enabled. Writes a per-descendant QC summary CSV.");

        forceTooltipSupport(outputField);
        forceTooltipSupport(siteWeightsField);
        forceTooltipSupport(apobecRateField);
        forceTooltipSupport(qcMetricBox);
        forceTooltipSupport(qcThresholdField);
        forceTooltipSupport(codonAwareBox);
        forceTooltipSupport(apobecBox);
        forceTooltipSupport(qcCheckBox);
        forceTooltipSupport(qcWriteCsvBox);

        for (java.awt.Component comp : form.getComponents()) {
            if (comp instanceof JLabel) {
                JLabel label = (JLabel) comp;
                String labelText = label.getText();
                if (labelText == null) {
                    continue;
                }
                if (labelText.startsWith("Output FASTA")) {
                    label.setToolTipText(outputField.getToolTipText());
                    forceTooltipSupport(label);
                } else if (labelText.startsWith("Site weights")) {
                    label.setToolTipText(siteWeightsField.getToolTipText());
                    forceTooltipSupport(label);
                } else if (labelText.startsWith("APOBEC rate")) {
                    label.setToolTipText(apobecRateField.getToolTipText());
                    forceTooltipSupport(label);
                } else if (labelText.startsWith("QC metric")) {
                    label.setToolTipText(qcMetricBox.getToolTipText());
                    forceTooltipSupport(label);
                } else if (labelText.startsWith("QC threshold")) {
                    label.setToolTipText(qcThresholdField.getToolTipText());
                    forceTooltipSupport(label);
                }
            }
        }
    }

    private void updateApobecFieldState() {
        apobecRateField.setEnabled(apobecBox.isSelected());
    }

    private void updateQcFieldState() {
        boolean qcEnabled = qcCheckBox.isSelected();
        qcMetricBox.setEnabled(qcEnabled);
        qcThresholdField.setEnabled(qcEnabled);
        qcWriteCsvBox.setEnabled(qcEnabled);
    }

    private void addRow(JPanel panel, GridBagConstraints c, int row,
                        String label, JComponent field, JComponent extra) {
        c.gridx = 0;
        c.gridy = row;
        c.weightx = 0;
        c.gridwidth = 1;
        panel.add(new JLabel(label), c);

        c.gridx = 1;
        c.weightx = 1.0;
        panel.add(field, c);

        if (extra != null) {
            c.gridx = 2;
            c.weightx = 0;
            panel.add(extra, c);
        }
    }

    private void updateDefaultOutputFromInput() {
        String inputPath = inputField.getText().trim();
        if (inputPath.isEmpty()) {
            return;
        }

        outputField.setText(defaultOutputFileForInput(new File(inputPath)).getAbsolutePath());
    }

    private File defaultOutputFileForInput(File inputFile) {
        String name = inputFile.getName();
        int dot = name.lastIndexOf('.');
        String base = (dot > 0) ? name.substring(0, dot) : name;

        File parent = inputFile.getParentFile();
        if (parent == null) {
            parent = new File(".");
        }

        return new File(parent, base + "Mut.fasta");
    }

    private void chooseFile(JTextField target, boolean saveDialog, boolean optionalTextFile) {
        JFileChooser chooser = (lastDirectory != null) ? new JFileChooser(lastDirectory) : new JFileChooser();

        int result = saveDialog ? chooser.showSaveDialog(this) : chooser.showOpenDialog(this);
        if (result == JFileChooser.APPROVE_OPTION) {
            File chosen = chooser.getSelectedFile();
            target.setText(chosen.getAbsolutePath());

            File parentDir = chosen.getParentFile();
            if (parentDir != null) {
                lastDirectory = parentDir;
            }

            if (!saveDialog && !optionalTextFile) {
                updateDefaultOutputFromInput();
            }
        }
    }

    private void cancelRun() {
        if (currentWorker != null && !currentWorker.isDone()) {
            currentWorker.cancel(true);
            System.out.println("Cancellation requested");
            logArea.append("Cancellation requested.\n");
            runButton.setEnabled(true);
            cancelButton.setEnabled(false);
        }
    }

    private void saveLog() {
        JFileChooser chooser = (lastDirectory != null) ? new JFileChooser(lastDirectory) : new JFileChooser();
        chooser.setDialogTitle("Save DUNES Log");

        int result = chooser.showSaveDialog(this);
        if (result != JFileChooser.APPROVE_OPTION) {
            return;
        }

        File selectedFile = chooser.getSelectedFile();
        File parentDir = selectedFile.getParentFile();
        if (parentDir != null) {
            lastDirectory = parentDir;
        }
        String path = selectedFile.getAbsolutePath();

        if (!path.toLowerCase().endsWith(".txt")) {
            selectedFile = new File(path + ".txt");
        }

        if (selectedFile.exists()) {
            int choice = JOptionPane.showConfirmDialog(
                    this,
                    "The log file already exists:\n" + selectedFile.getAbsolutePath() + "\n\nOverwrite it?",
                    "Confirm overwrite",
                    JOptionPane.YES_NO_OPTION,
                    JOptionPane.WARNING_MESSAGE
            );

            if (choice != JOptionPane.YES_OPTION) {
                return;
            }
        }

        try (FileWriter writer = new FileWriter(selectedFile)) {
            writer.write(logArea.getText());
            JOptionPane.showMessageDialog(
                    this,
                    "Log saved to:\n" + selectedFile.getAbsolutePath(),
                    "Log saved",
                    JOptionPane.INFORMATION_MESSAGE
            );
        } catch (Exception ex) {
            JOptionPane.showMessageDialog(
                    this,
                    "Could not save log:\n" + ex.getMessage(),
                    "Save error",
                    JOptionPane.ERROR_MESSAGE
            );
        }
    }

    private void clearLog() {
        if (!logArea.getText().trim().isEmpty()) {
            int choice = JOptionPane.showConfirmDialog(
                    this,
                    "Clear the current log?",
                    "Confirm clear log",
                    JOptionPane.YES_NO_OPTION,
                    JOptionPane.QUESTION_MESSAGE
            );

            if (choice != JOptionPane.YES_OPTION) {
                return;
            }
        }

        logArea.setText("");
        progressBar.setValue(0);
        progressBar.setString("0%");
    }

    private void resetForm() {
        int choice = JOptionPane.showConfirmDialog(
                this,
                "Reset all form fields to the default DUNES settings?",
                "Confirm reset form",
                JOptionPane.YES_NO_OPTION,
                JOptionPane.QUESTION_MESSAGE
        );

        if (choice != JOptionPane.YES_OPTION) {
            return;
        }

        inputField.setText("");
        outputField.setText("");
        siteWeightsField.setText("");

        mutationField.setText("0.004");
        yearsField.setText("2");
        mutantsField.setText("10");
        modelBox.setSelectedItem("HIV");
        codonAwareBox.setSelected(false);
        apobecBox.setSelected(false);
        apobecRateField.setText("0.02");
        qcCheckBox.setSelected(false);
        qcMetricBox.setSelectedItem("normalized_nt_distance");
        qcThresholdField.setText("0.01");
        qcWriteCsvBox.setSelected(true);

        updateApobecFieldState();
        updateQcFieldState();

        progressBar.setValue(0);
        progressBar.setString("0%");

        logArea.append("Form reset to default settings.\n");
    }

    private void showHelpDialog() {
        String helpText =
                "DUNES GUI\n\n" +
                "Purpose:\n" +
                "Generate synthetic descendant nucleotide sequences from an input FASTA file.\n\n" +

                "Main settings:\n" +
                "- Input FASTA: source sequences to mutate.\n" +
                "- Output FASTA (optional): file where descendants will be written.\n" +
                "  If left blank, DUNES writes <input base>Mut.fasta next to the input FASTA.\n" +
                "- Mutation rate: baseline substitutions per nucleotide per year.\n" +
                "- Years: evolutionary time used to calculate mutation probability.\n" +
                "- Mutants / sequence: number of descendants generated for each parent sequence.\n" +
                "- Model:\n" +
                "    simple = approximately even non-self substitutions\n" +
                "    HIV = HIV-like biased substitution model\n\n" +

                "Options:\n" +
                "- Codon-aware: filters mutations using a simple coding-sequence rule set.\n" +
                "  Useful for testing, but usually not needed for subtype benchmarking.\n" +
                "- APOBEC: applies stronger G->A bias to a subset of descendants.\n" +
                "- APOBEC rate (optional): probability that a descendant uses the APOBEC-enhanced mutation model.\n" +
                "- Site weights (optional): text file with one numeric value per line.\n" +
                "  The number of lines must match the sequence length.\n" +
                "- Enable QC distance check: compare each descendant to its parent after generation.\n" +
                "- QC metric (optional):\n" +
                "    normalized_nt_distance = differences / sequence length\n" +
                "    raw_nt_differences = total differing nucleotide positions\n" +
                "- QC threshold (optional): maximum allowed distance for pass/fail reporting.\n" +
                "- Write QC CSV: writes a per-descendant QC summary file.\n" +
                "- Site-weights validation is automatic; there is no separate matching/mismatching option.\n\n" +

                "Recommended HIV subtype benchmarking defaults:\n" +
                "- mutation rate = 0.004\n" +
                "- years = 2\n" +
                "- mutants / sequence = 10\n" +
                "- model = HIV\n" +
                "- codon-aware = off\n" +
                "- APOBEC = off\n" +
                "- site weights = blank\n" +
                "- QC disabled by default; if enabled, start with normalized distance threshold = 0.01\n\n" +

                "Notes:\n" +
                "- Ambiguity-containing sequences are supported.\n" +
                "- Mutated ambiguity sites resolve to canonical bases when a mutation occurs.\n" +
                "- Large runs may take time; use Cancel if needed.\n" +
                "- Save Log can export run details to a text file.\n";

        JTextArea helpArea = new JTextArea(helpText, 24, 70);
        helpArea.setWrapStyleWord(true);
        helpArea.setLineWrap(true);
        helpArea.setEditable(false);
        helpArea.setCaretPosition(0);

        JScrollPane scrollPane = new JScrollPane(helpArea);

        JOptionPane.showMessageDialog(
                this,
                scrollPane,
                "DUNES Help / About",
                JOptionPane.INFORMATION_MESSAGE
        );
    }

    private void openManual() {
        try {
            File appDir = getApplicationDirectory();

            File[] candidates = new File[] {
                    new File(appDir, "DUNES_User_Manual.docx"),
                    new File(appDir, "DUNES User Manual.docx"),
                    new File(appDir, "DUNES_User_Manual.pdf"),
                    new File(appDir, "DUNES User Manual.pdf")
            };

            File manualFile = null;
            for (File f : candidates) {
                if (f.exists()) {
                    manualFile = f;
                    break;
                }
            }

            if (manualFile == null) {
                JOptionPane.showMessageDialog(
                        this,
                        "Could not find the user manual next to the DUNES jar.\n\nChecked folder:\n"
                                + appDir.getAbsolutePath(),
                        "Manual not found",
                        JOptionPane.WARNING_MESSAGE
                );
                return;
            }

            if (!Desktop.isDesktopSupported()) {
                JOptionPane.showMessageDialog(
                        this,
                        "Desktop file opening is not supported on this system.\n\nManual path:\n" + manualFile.getAbsolutePath(),
                        "Open manual unavailable",
                        JOptionPane.WARNING_MESSAGE
                );
                return;
            }

            Desktop.getDesktop().open(manualFile);

        } catch (Exception ex) {
            ex.printStackTrace();
            JOptionPane.showMessageDialog(
                    this,
                    "Could not open the manual:\n" + ex.toString(),
                    "Open manual error",
                    JOptionPane.ERROR_MESSAGE
            );
        }
    }

    private void runDunes() {
        final File inFile;
        final File outFile;
        final double mutationRate;
        final double years;
        final int n;
        final String model;
        final boolean codonAware;
        final boolean apobec;
        final double apobecRate;
        final File siteWeightsFile;
        final boolean qcEnabled;
        final String qcMetric;
        final double qcThreshold;
        final boolean qcWriteCsv;

        try {
            String inputPath = inputField.getText().trim();
            String outputPath = outputField.getText().trim();
            String weightsPath = siteWeightsField.getText().trim();

            if (inputPath.isEmpty()) {
                throw new IllegalArgumentException("Please select an input FASTA file.");
            }
            inFile = new File(inputPath);
            if (outputPath.isEmpty()) {
                outFile = defaultOutputFileForInput(inFile);
                outputField.setText(outFile.getAbsolutePath());
            } else {
                outFile = new File(outputPath);
            }

            if (!inFile.exists()) {
                throw new IllegalArgumentException("Input FASTA file does not exist.");
            }

            if (outFile.exists()) {
                int choice = JOptionPane.showConfirmDialog(
                        this,
                        "The output file already exists:\n" + outFile.getAbsolutePath() + "\n\nOverwrite it?",
                        "Confirm overwrite",
                        JOptionPane.YES_NO_OPTION,
                        JOptionPane.WARNING_MESSAGE
                );

                if (choice != JOptionPane.YES_OPTION) {
                    return;
                }
            }

            if (inFile.length() > MAX_SUGGESTED_FILE_SIZE) {
                int choice = JOptionPane.showConfirmDialog(
                        this,
                        "The input FASTA is larger than 100 MB.\n"
                                + "Large files may take a long time or use a lot of memory.\n\n"
                                + "Do you want to continue?",
                        "Large input warning",
                        JOptionPane.YES_NO_OPTION,
                        JOptionPane.WARNING_MESSAGE
                );

                if (choice != JOptionPane.YES_OPTION) {
                    return;
                }
            }

            mutationRate = Double.parseDouble(mutationField.getText().trim());
            years = Double.parseDouble(yearsField.getText().trim());
            n = Integer.parseInt(mutantsField.getText().trim());
            String selectedModel = (String) modelBox.getSelectedItem();
            model = "HIV".equalsIgnoreCase(selectedModel) ? "hiv" : selectedModel.toLowerCase();
            codonAware = codonAwareBox.isSelected();
            apobec = apobecBox.isSelected();
            apobecRate = apobec ? Double.parseDouble(apobecRateField.getText().trim()) : 0.0;

            qcEnabled = qcCheckBox.isSelected();
            qcMetric = (String) qcMetricBox.getSelectedItem();
            qcWriteCsv = qcWriteCsvBox.isSelected();
            qcThreshold = qcEnabled ? Double.parseDouble(qcThresholdField.getText().trim()) : 0.0;

            if (mutationRate < 0.0) {
                throw new IllegalArgumentException("Mutation rate must be >= 0.");
            }
            if (years < 0.0) {
                throw new IllegalArgumentException("Years must be >= 0.");
            }
            if (n < 1) {
                throw new IllegalArgumentException("Mutants / sequence must be at least 1.");
            }
            if (apobecRate < 0.0 || apobecRate > 1.0) {
                throw new IllegalArgumentException("APOBEC rate must be between 0 and 1.");
            }
            if (qcEnabled && qcThreshold < 0.0) {
                throw new IllegalArgumentException("QC threshold must be >= 0.");
            }

            if (weightsPath.isEmpty()) {
                siteWeightsFile = null;
            } else {
                siteWeightsFile = new File(weightsPath);
                if (!siteWeightsFile.exists()) {
                    throw new IllegalArgumentException("Site-weights file does not exist.");
                }

                if (siteWeightsFile.length() > MAX_SUGGESTED_FILE_SIZE) {
                    int choice = JOptionPane.showConfirmDialog(
                            this,
                            "The site-weights file is larger than 100 MB.\n"
                                    + "This may take a long time to load.\n\n"
                                    + "Do you want to continue?",
                            "Large file warning",
                            JOptionPane.YES_NO_OPTION,
                            JOptionPane.WARNING_MESSAGE
                    );

                    if (choice != JOptionPane.YES_OPTION) {
                        return;
                    }
                }
            }

        } catch (Exception ex) {
            JOptionPane.showMessageDialog(this, ex.getMessage(), "Input error", JOptionPane.ERROR_MESSAGE);
            return;
        }

        System.out.println("Run requested");
        logArea.append("DEBUG: Run requested\n");

        System.out.println("Creating background worker");
        logArea.append("DEBUG: Creating background worker\n");

        runButton.setEnabled(false);
        cancelButton.setEnabled(true);
        progressBar.setValue(0);
        progressBar.setString("0%");
        logArea.append("Running DUNES at " + LocalDateTime.now().format(TS_FORMAT) + "...\n");
        logArea.append("Input: " + inFile.getAbsolutePath() + "\n");
        logArea.append("Output: " + outFile.getAbsolutePath() + "\n");
        logArea.append("Model: " + ("hiv".equalsIgnoreCase(model) ? "HIV" : model) + "\n");
        logArea.append("Mutation rate: " + mutationRate + "\n");
        logArea.append("Years: " + years + "\n");
        logArea.append("Mutants / sequence: " + n + "\n");
        logArea.append("Codon-aware: " + codonAware + "\n");
        logArea.append("APOBEC: " + apobec + "\n");
        if (apobec) {
            logArea.append("APOBEC rate: " + apobecRate + "\n");
        }
        if (siteWeightsFile != null) {
            logArea.append("Site weights: " + siteWeightsFile.getAbsolutePath() + "\n");
        }
        logArea.append("QC enabled: " + qcEnabled + "\n");
        if (qcEnabled) {
            logArea.append("QC metric: " + qcMetric + "\n");
            logArea.append("QC threshold: " + qcThreshold + "\n");
            logArea.append("Write QC CSV: " + qcWriteCsv + "\n");
        }
        logArea.append("\n");

        final DunesEngine engine = new DunesEngine();

        currentWorker = new SwingWorker<Void, Void>() {
            @Override
            protected Void doInBackground() throws Exception {
                System.out.println("Background worker started");
                SwingUtilities.invokeLater(() -> logArea.append("DEBUG: Background worker started\n"));

                if (isCancelled()) {
                    return null;
                }

                System.out.println("Calling DunesEngine.run()");
                SwingUtilities.invokeLater(() -> logArea.append("DEBUG: Calling DunesEngine.run()\n"));

                engine.run(
                        inFile,
                        outFile,
                        mutationRate,
                        years,
                        n,
                        model,
                        codonAware,
                        apobec,
                        apobecRate,
                        siteWeightsFile,
                        (message, currentParent, totalParents) -> SwingUtilities.invokeLater(() -> {
                            logArea.append(message + "\n");
                            int percent = (totalParents > 0) ? (int) Math.round((currentParent * 100.0) / totalParents) : 0;
                            progressBar.setValue(percent);
                            progressBar.setString(percent + "%");
                        }),
                        qcEnabled,
                        qcMetric,
                        qcThreshold
                );
                return null;
            }

            @Override
            protected void done() {
                runButton.setEnabled(true);
                cancelButton.setEnabled(false);

                try {
                    if (isCancelled()) {
                        progressBar.setString("Cancelled");
                        logArea.append("Run cancelled.\n\n");
                        return;
                    }

                    get();

                    logArea.append("Done. Output written to: " + outFile.getAbsolutePath() + "\n");
                    logArea.append("Summary: processed " + engine.getLastParentCount()
                            + " parent sequence(s) and wrote "
                            + engine.getLastDescendantCount() + " descendant sequence(s).\n");
                    logArea.append("Run settings: model=" + ("hiv".equalsIgnoreCase(model) ? "HIV" : model)
                            + ", mutation_rate=" + mutationRate
                            + ", years=" + years
                            + ", mutants_per_sequence=" + n
                            + ", codon_aware=" + codonAware
                            + ", apobec=" + apobec
                            + ", apobec_rate=" + apobecRate
                            + ", site_weights=" + (siteWeightsFile == null ? "none" : siteWeightsFile.getAbsolutePath())
                            + ", qc_enabled=" + qcEnabled
                            + ", qc_metric=" + qcMetric
                            + ", qc_threshold=" + qcThreshold
                            + "\n");

                    if (qcEnabled) {
                        logArea.append("QC summary: " + engine.getLastQcPassCount()
                                + " pass, " + engine.getLastQcFailCount() + " fail.\n");

                        if (qcWriteCsv) {
                            File qcCsv = DunesEngine.defaultQcCsvPath(outFile);
                            engine.writeLastQcResultsCsv(qcCsv);
                            logArea.append("QC CSV written to: " + qcCsv.getAbsolutePath() + "\n");
                        }
                    }

                    logArea.append("Completed at " + LocalDateTime.now().format(TS_FORMAT) + "\n\n");

                    progressBar.setValue(100);
                    progressBar.setString("Complete");
                } catch (Exception ex) {
                    Throwable root = ex;
                    if (ex instanceof java.util.concurrent.ExecutionException && ex.getCause() != null) {
                        root = ex.getCause();
                    }

                    root.printStackTrace();

                    progressBar.setString("Error");
                    logArea.append("ERROR: " + root.toString() + "\n");

                    if (root.getMessage() != null && !root.getMessage().trim().isEmpty()) {
                        logArea.append("DETAIL: " + root.getMessage() + "\n");
                    }

                    logArea.append("\n");

                    JOptionPane.showMessageDialog(
                            DunesGui.this,
                            "Run failed: " + root.toString(),
                            "Error",
                            JOptionPane.ERROR_MESSAGE
                    );
                } finally {
                    currentWorker = null;
                }
            }
        };

        currentWorker.execute();
    }
}

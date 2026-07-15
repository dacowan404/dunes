# DUNES

DUNES introduces random mutations into DNA FASTA sequences and writes synthetic descendant sequences.

 - Download full version with GUI: [dunes-full.jar](https://github.com/dacowan404/dunes/releases/download/0.2/dunes-full.jar)
 - Download command line only version: [dunes-cli.jar](https://github.com/dacowan404/dunes/releases/download/0.2/dunes-cli.jar)
 - Download user manual: [Dunes User Manual](https://github.com/dacowan404/dunes/releases/download/0.2/DUNES_User_Manual.docx)

## Requirements

- Java 11 or newer is recommended.
- The GUI shows a warning when it is started with a Java runtime older than 11.
- Maven is required only when compiling from source.

## Build

```bash
mvn clean package
```

The package step creates two executable jars with dependencies:

- `target/dunes-full.jar`: full application with GUI and command-line support.
- `target/dunes-cli.jar`: command-line-only jar without the GUI classes.

Maven also creates the thin project jar, but the two jars above are the distributable outputs because they include dependencies.

## Run The GUI

```bash
java -jar target/dunes-full.jar
```

Running `dunes-full.jar` with no arguments opens the GUI. Optional GUI fields are labeled as optional:

- Output FASTA: defaults to `<input base>Mut.fasta` next to the input FASTA when left blank.
- Site weights: leave blank to use equal weight at every site.
- APOBEC rate: used only when APOBEC is enabled.
- QC metric and QC threshold: used only when QC is enabled.

## Run From The Command Line

```bash
java -jar target/dunes-cli.jar -i input.fasta
```

The full jar can also run command-line jobs when arguments are supplied:

```bash
java -jar target/dunes-full.jar -i input.fasta --model HIV --mutants-number 10
```

## Command-Line Options

```text
-h, --help
    Show help and exit.

-V, --version
    Print version information and exit.

-i, --inFile=FILE
    Required input FASTA file.

-o, --outFile=FILE
    Optional output FASTA file. Defaults to <input base>Mut.fasta.

-m, --mutation-rate=VALUE
    Mutation rate per nucleotide per year. Default: 0.004.

-y, --years=VALUE
    Years of evolution. Default: 2.0.

-n, --mutants-number=VALUE
    Number of descendants to generate per input sequence. Default: 10.

--model=VALUE
    Evolution model: simple or HIV. Default: hiv.

--codon-aware
    Enable the codon-aware acceptance filter.

--apobec
    Enable APOBEC-like G-to-A hypermutation for a subset of descendants.

--apobec-rate=VALUE
    Probability that a descendant is APOBEC affected. Default: 0.02.

--site-weights=FILE
    Optional text file with one numeric weight per line. The number of lines must match the sequence length.

--qc-enable
    Enable QC distance checks between each parent and descendant.

--qc-metric=VALUE
    QC metric: normalized_nt_distance or raw_nt_differences. Default: normalized_nt_distance.

--qc-threshold=VALUE
    Maximum allowed QC distance for pass/fail reporting. Default: 0.01.

--qc-write-csv
    Write a QC CSV summary next to the output FASTA. Enabled by default when QC results are written.
```

## Codon-Aware And APOBEC Behavior

- Codon-aware mode evaluates proposed substitutions in codons starting at the first nucleotide. Synonymous changes are accepted, changes that introduce a stop codon are rejected, conservative amino-acid changes have a 50% chance of being accepted, and other amino-acid changes have a 20% chance of being accepted. Codons containing ambiguity or gap characters, and incomplete codons at the end of a sequence, are not filtered.
- APOBEC mode uses an enhanced G-to-A mutation bias only with the `hiv` model. Enabling APOBEC with the `simple` model does not change the simple mutation matrix. For the HIV model, `--apobec-rate` controls the probability that each descendant uses the APOBEC-enhanced matrix.

## Examples

Generate 10 HIV-model descendants per input sequence with default settings:

```bash
java -jar target/dunes-cli.jar -i input.fasta
```

Generate simple-model descendants and write to a specific output file:

```bash
java -jar target/dunes-cli.jar -i input.fasta -o output.fasta --model simple -n 25
```

Run with APOBEC and QC enabled:

```bash
java -jar target/dunes-cli.jar -i input.fasta --apobec --apobec-rate 0.05 --qc-enable --qc-threshold 0.01
```

## Mutation Models

Rows indicate the original base and columns indicate the mutated base. These probabilities are applied only after a mutation event has been selected for a nucleotide.

### `simple`

| From \ To | A | C | G | T |
| --- | ---: | ---: | ---: | ---: |
| A | 0.0 | 0.33 | 0.33 | 0.33 |
| C | 0.33 | 0.0 | 0.33 | 0.33 |
| G | 0.33 | 0.33 | 0.0 | 0.33 |
| T | 0.33 | 0.33 | 0.33 | 0.0 |

### `hiv`

| From \ To | A | C | G | T |
| --- | ---: | ---: | ---: | ---: |
| A | 0.0 | 0.15 | 0.70 | 0.15 |
| C | 0.15 | 0.0 | 0.15 | 0.70 |
| G | 0.80 | 0.10 | 0.0 | 0.10 |
| T | 0.15 | 0.70 | 0.15 | 0.0 |

![DUNES](https://github.com/Sergey-Knyazev/dunes/blob/master/DUNES%20icon_2.svg)

# DUNES
A command line tool for introducing random mutations into DNA sequences.

## Download or Compilation

- Download:
[dunes.jar](https://github.com/dacowan404/dunes/releases/tag/0.1.2/dunes-jar-with-dependencies.jar)

- Compilation:
```bash
mvn clean install
```

## Running

```bash
java -jar dunes.jar
```

## Help

```bash
java -jar dunes.jar -h

Usage: dunes [-hV] -i=FILE [-m=m] [-n=n] [-o=<outputFolder>] [-y=y] [-d=d]
  -h, --help               Show this help message and exit.
  -d, --distribution <distribution>       mutation model to use: 'simple' (equal-probability) or 'hiv' (HIV-like) (default: 'simple')
  -i, --inFile=FILE        input fasta file with viral sequences to be mutated
  -m, --mutation-rate=m    mutation rate in substitutions per nucleotide per year
                             (s/n/y) (default: 4.1e-3)
  -n, --mutants-number=n   number of mutants for a strain (default: 1)
  -o, --outFile=<outputFolder>
                           fasta file with mutated sequences (default: "<inFile
                             baseMut.fasta>")
  -V, --version            Print version information and exit.
  -y, --years=y            years of evolution (default: 1)
```

## Distribution Matrices

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

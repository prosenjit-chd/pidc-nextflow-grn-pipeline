# 🔬 PIDC Gene Regulatory Network Inference Pipeline (Nextflow)

This project implements the **PIDC (Partial Information Decomposition and Context)** algorithm for gene regulatory network (GRN) inference from single-cell RNA-seq data using a **Nextflow DSL2 pipeline**. It supports `.h5ad` format as input and outputs a ranked list of gene-gene regulatory relationships.

The pipeline is based on [NetworkInference.jl](https://github.com/Tchanders/NetworkInference.jl) and follows best practices from the [Beeline GRN Benchmarking Suite](https://github.com/Murali-group/Beeline).

---

## 📁 Project Structure

```

pidc\_nextflow/
├── main.nf                        # Main Nextflow pipeline script
├── nextflow\.config                # Configuration file (parameters + conda)
├── crosscheck.py                  # Gene overlap checker + visualizer
├── results/                       # Output folder (auto-created)
├── modules/
│   └── grn/
│       └── pidc/
│           ├── environment.yml   # Conda environment file for scanpy, pandas
│           ├── runPIDC.jl        # Julia script (PIDC implementation)
│           └── run\_pidc.sh       # Shell script to call Julia
└── filtered\_placeholder.h5ad     # Example synthetic input data

```

---

## 🔧 Setup Instructions

### 🧬 Prerequisites

Make sure you have the following installed:

- [Nextflow](https://www.nextflow.io/) ≥ 22.x
- [Miniconda/Conda](https://docs.conda.io/en/latest/)
- [Julia](https://julialang.org/) ≥ 1.6
- Python 3 with:
  - `scanpy`
  - `pandas`
  - `matplotlib`

### ✅ Environment Setup

```bash
git clone https://github.com/YOUR_USERNAME/pidc_nextflow.git
cd pidc_nextflow

# Create Conda environment
conda env create -f modules/grn/pidc/environment.yml -p ./test_env

# Activate it
conda activate ./test_env
```

---

## 🚀 Run the PIDC Pipeline

```bash
NXF_CONDA_CACHEDIR=./.conda nextflow run main.nf \
  --expression_h5ad filtered_placeholder.h5ad \
  --prefix pidc_output_final \
  -profile conda
```

**Result:**
It will create a CSV file in the `results/` folder:

```
results/pidc_output_final_rankedEdges.csv
```

---

## 📤 Output Format

PIDC outputs a **ranked edge list** in the following format:

```
Gene1    Gene2    Score
```

**Example:**

```
DENND2D    RCN2     1.9742
RCN2       DENND2D  1.9742
PRAF2      DENND2D  1.9728
```

> 🔁 The output includes **bidirectional pairs**.

---

## 📊 Post-processing Visualization

Use the `crosscheck.py` script to compare input vs output genes:

```bash
python3 crosscheck.py
```

This will generate:

```
gene_comparison.png
```

### 🧾 Sample Output Categories

- ✅ **Only in Input**: Genes in `.h5ad` but not in output
- ✅ **Only in Output**: Genes in output but not in `.h5ad`
- ✅ **Common Genes**: Successfully inferred genes

---

## 📚 Reference

- **Chan, Stumpf, Babtie (2017)**
  _Gene Regulatory Network Inference from Single-Cell Data Using Multivariate Information Measures_
  [Cell Systems, Vol 5, Issue 3](https://doi.org/10.1016/j.cels.2017.08.002)

- [NetworkInference.jl](https://github.com/Tchanders/NetworkInference.jl)

- [Beeline Benchmark Suite](https://github.com/Murali-group/Beeline)

---

## 👨‍💻 Author

**Prosenjit Chowdhury**
M.Sc. Artificial Intelligence – FAU Erlangen-Nürnberg
Working Student @ SAP ERP PCX
🌍 Erlangen, Germany
🔗 GitHub: [@prosenjit-chowdhury](https://github.com/prosenjit-chowdhury)

---

## 🧠 License

This project is licensed under the MIT License. See [`LICENSE`](LICENSE) for details.

```

```

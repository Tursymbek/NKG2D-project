# Taiji 1.3 — Complete Setup and Run Guide
## NKG2D Project: ATAC-seq + RNA-seq Transcription Factor Analysis

---

## System Specifications

| Component | Details |
|-----------|---------|
| Device | Legion 82B5 |
| CPU | AMD Ryzen 5 4600H with Radeon Graphics (6 cores / 12 threads, 3.00 GHz) |
| RAM | 32.0 GB |
| GPU | NVIDIA GeForce GTX 1650 Ti (4 GB) |
| Storage | 943 GB total |
| OS | Windows 11 (26200.8457), 64-bit |
| WSL Version | WSL 2.6.1.0 |
| WSL Kernel | 6.6.87.2-1 |

---

## Software Versions

| Package | Version |
|---------|---------|
| OS (WSL) | Ubuntu 24.04.3 LTS (Noble) |
| WSL Kernel | 6.6.87.2-microsoft-standard-WSL2 |
| Taiji | v1.3.0 (taiji-Ubuntu-x86_64) |
| bedtools | v2.31.1 |
| samtools | v1.19.2 |
| sqlite3 | 3.51.0 (2025-11-04) |
| Python | 3.13.9 |
| wget | 1.21.4 |
| curl | 8.5.0 |

---

## WSL2 Configuration

File location: `C:\Users\Lenovo\.wslconfig`

```ini
[wsl2]
processors=12
memory=28GB
swap=4GB
localhostForwarding=true
```

Apply changes:
```powershell
wsl --shutdown
```

---

## Project Structure

```
/mnt/d/NKG2D project/TAIJIv1.3/
├── ATAC/                        # only IDR thresholded peaks are used
│   ├── ENCFF018OJP.bed.gz       # A549 (NarrowPeak)
│   ├── ENCFF020EUE.bed.gz
│   └── ... (81 files total)
├── RNA/
│   ├── A549_gene_quant.tsv      # gene + TPM, no header
│   ├── HCT_116_gene_quant.tsv
│   ├── PC_3_gene_quant.tsv
│   ├── Hep_G2_gene_quant.tsv
│   └── NCI_H929_gene_quant.tsv
├── reference/
│   └── hg38.fa
├── gencode.v49.basic.annotation.gtf
├── config.yaml
├── samples.yml
├── sciflow.db                   # Taiji cache — DO NOT DELETE unless necessary
└── results/
```

---

## Sample Groups

| Cell Line | Group | Data Type | Replicates |
|-----------|-------|-----------|------------|
| A549 | High | ATAC-seq | 36 |
| HCT116 | High | ATAC-seq | 27 |
| PC-3 | High | ATAC-seq | 3 |
| HepG2 | Low | ATAC-seq | 9 |
| NCI-H929 | Low | ATAC-seq | 6 |
| A549 | High | RNA-seq | 1 (GeneQuant) |
| HCT116 | High | RNA-seq | 1 (GeneQuant) |
| PC-3 | High | RNA-seq | 1 (GeneQuant) |
| HepG2 | Low | RNA-seq | 1 (GeneQuant) |
| NCI-H929 | Low | RNA-seq | 1 (GeneQuant) |

---

## config.yaml

```yaml
assembly: hg38

genome: "/mnt/d/NKG2D project/TAIJIv1.3/reference/hg38.fa"

annotation: "/mnt/d/NKG2D project/TAIJIv1.3/gencode.v49.basic.annotation.gtf"

input: samples.yml

output_dir: results

threads: 8

network:
  enhancerPromoterDistance: 500000
```

---

## samples.yml (format)

Key rules:
- ATAC-seq files: use `format: NarrowPeak` (NOT in tags)
- RNA-seq GeneQuant files: use `tags: ['GeneQuant']` (NOT in format)
- No empty lines
- Indentation with spaces only (no tabs)

```yaml
ATAC-seq:
  - id: A549_ATAC
    group: High
    replicates:
      - rep: 1
        files:
          - path: ATAC/ENCFF018OJP.bed.gz
            format: NarrowPeak

RNA-seq:
  - id: A549_RNA
    group: High
    replicates:
      - rep: 1
        files:
          - path: RNA/A549_gene_quant.tsv
            tags: ['GeneQuant']
```

---

## Installation

```bash
# 1. Install dependencies
sudo apt-get update
sudo apt-get install -y bedtools samtools sqlite3 wget curl

# 2. Download Taiji binary
curl -L https://github.com/Taiji-pipeline/Taiji/releases/latest/download/taiji-Ubuntu-x86_64 \
  -o ~/taiji-Ubuntu-x86_64
chmod +x ~/taiji-Ubuntu-x86_64

# 3. Verify
~/taiji-Ubuntu-x86_64 --help
```

---

## Running Taiji

### Navigate to project folder
```bash
cd /mnt/d/NKG2D\ project/TAIJIv1.3/
```

### Step 1 — 7: ATAC_Find_TFBS (CPU-intensive)
Use maximum threads:
```bash
~/taiji-Ubuntu-x86_64 run --config config.yaml -n 12 +RTS -N12
```

### Step 8: Create_Linkage and beyond (memory-intensive)
Reduce threads to avoid OOM crash:
```bash
~/taiji-Ubuntu-x86_64 run --config config.yaml -n 1 +RTS -N1
```

> **Note:** Taiji resumes from where it stopped — it is safe to stop and restart at any time.

---

## Pipeline Steps and Expected Duration

| Step | Description | 
|------|-------------|
| Read_Input | Parse samples.yml |
| Download_Data | Copy/verify local files |
| Make_Index | Build genome index |
| ATAC_Merge_Peaks | Merge peaks across replicates |
| **ATAC_Find_TFBS** | **Scan TF motifs in open chromatin** |
| RNA_Make_Expr_Table | Build expression matrix |
| Create_Linkage | Link ATAC peaks to gene promoters |
| Compute_Ranks | PageRank for each TF |
| Output_Ranks | Write GeneRanks.tsv |

---

## Output Files

```
results/
├── GeneRanks.tsv              # Main result: TF PageRank per group
├── GeneRanks_PValues.tsv      # Statistical significance
├── Network/
│   ├── High/
│   │   ├── edges_combined.csv # TF → target gene edges
│   │   ├── edges_binding.csv  # ATAC binding sites
│   │   └── nodes.csv          # Gene expression weights
│   └── Low/
│       ├── edges_combined.csv
│       ├── edges_binding.csv
│       └── nodes.csv
├── ATACSeq/
│   ├── openChromatin.bed.gz   # Merged open chromatin regions
│   └── TFBS/                  # TF binding site files per sample
└── RNASeq/
    └── expression_profile.tsv # Averaged expression per group
```

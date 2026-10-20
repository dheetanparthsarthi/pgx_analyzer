---
title: PGx Population Pharmacogenomics Platform
emoji: 🧬
colorFrom: blue
colorTo: indigo
sdk: streamlit
sdk_version: 1.41.1
app_file: app.py
pinned: false
---

# 🧬 Automated Population Pharmacogenomics Analysis Platform

A population pharmacogenomics platform for genotype and allele frequencies, Hardy–Weinberg Equilibrium (HWE), regional/gender stratification, and rule-configured diplotype/phenotype classification. Marker-level analysis is data-driven; clinical translation is available for genes with validated rule bundles.

---

## 🎯 1. Platform Aim & What the Tool Does Exactly

Currently, population pharmacogenomic analysis involves studying the distribution of pharmacogenetic variants across diverse geographic and demographic cohorts. Manual spreadsheet calculations for genotype counts, allele frequencies, Hardy–Weinberg Equilibrium (HWE), demographic subgrouping, and diplotype assignment are time-consuming, prone to human error, and unscalable for growing patient cohorts.

**This platform automates the complete end-to-end analytical workflow:**
- **Automated Data Validation & QC**: Parses `.xlsx`, `.xls`, and `.csv` files, standardizing raw column headers, detecting duplicate Sample IDs, and quarantining invalid allele calls.
- **Genetic & Statistical Engine**: Instantly computes genotype frequencies, allele frequencies ($p$ and $q$), observed vs. expected counts, Chi-Square goodness-of-fit statistics ($\chi^2$), P-values, and Haldane exact test fallbacks.
- **Demographic & Regional Stratification**: Categorizes samples across 5 Indian geographic regions (**South India, North India, East India, West India, Central India**), 27+ individual state profiles, and Male vs. Female subgroups.
- **Rule-Based Pharmacogenomic Classifier**: Uses shared translation code with per-gene marker, diplotype, and phenotype rules. CYP2C19 and CYP2C9 are configured today; other genes receive population marker statistics until validated translation rules are supplied.
- **Multi-Format Exporters**: Generates publication-ready multi-tab Excel (`.xlsx`), summary CSV (`.csv`), and formal PDF (`.pdf`) summary reports with one click.

---

## ⚡ 2. Automated End-to-End Pipeline Architecture

The platform executes an 11-stage automated pipeline upon dataset upload:

```
  ┌────────────────────────┐
  │ Excel / CSV File Upload│ (.xlsx, .xls, .csv)
  └───────────┬────────────┘
              │
              ▼
  ┌────────────────────────┐
  │ Data Validation & QC   │ (Header normalization, duplicate detection, invalid allele quarantine)
  └───────────┬────────────┘
              │
              ▼
  ┌────────────────────────┐
  │ Genotype Normalization │ (Genotype string cleaning: G/A -> GA, T/C -> CT)
  └───────────┬────────────┘
              │
              ▼
  ┌────────────────────────┐
  │ Genotype Frequencies   │ (N_AA, N_Aa, N_aa counts, frequencies, percentages)
  └───────────┬────────────┘
              │
              ▼
  ┌────────────────────────┐
  │ Allele Frequencies     │ (Ref/Var allele counts, p and q frequencies)
  └───────────┬────────────┘
              │
              ▼
  ┌────────────────────────┐
  │ Hardy-Weinberg Engine  │ (Observed vs Expected, Chi-Square, P-value, Haldane Exact Test)
  └───────────┬────────────┘
              │
              ▼
  ┌────────────────────────┐
  │ Population Grouping    │ (South, North, East, West, Central India & Gender Subgroups)
  └───────────┬────────────┘
              │
              ▼
  ┌────────────────────────┐
  │ Star Allele Calling    │ (CYP2C19 *2, *3, *17 variant combination mapping)
  └───────────┬────────────┘
              │
              ▼
  ┌────────────────────────┐
  │ Diplotype Assignment   │ (*1/*1, *1/*2, *1/*17, *2/*17, *2/*2, *17/*17)
  └───────────┬────────────┘
              │
              ▼
  ┌────────────────────────┐
  │ Phenotype Calling      │ (CPIC classification: UM, RM, NM, IM, PM)
  └───────────┬────────────┘
              │
              ▼
  ┌────────────────────────┐
  │ Interactive Dashboard  │ (9-Tab Streamlit UI, Plotly charts, HTML warning badges)
  └───────────┬────────────┘
              │
              ▼
  ┌────────────────────────┐
  │ Multi-Format Exporters │ (Excel .xlsx, CSV .csv, PDF .pdf reports)
  └────────────────────────┘
```

---

## 💾 3. Database Architecture & Automated Batch Updates

### 3.1 Why SQLite Was Chosen
The platform utilizes **SQLite** (`database/db_manager.py` $\rightarrow$ `pharmacogenomics.db`) as its core persistence backend.
- **Zero Configuration Overhead**: Lightweight, serverless embedded database requiring zero database administration or external daemon processes.
- **ACID Compliance & Reliability**: Ensures atomic transactions, preventing database corruption during batch sample writes.
- **Local Zero-Latency Access**: Disk-backed persistence allows instant query execution and in-memory caching.
- **High Performance Scaling**: Optimized to handle initial cohorts of **1,000 samples**, scaling seamlessly to **1,200**, **2,000**, and **10,000+ samples**.

### 3.2 How Data Updates & Recalculate Automatically
1. **Deduplication Engine**: Every sample is assigned a primary key (`Sample_ID`). When a new sample batch is ingested, the engine verifies `Sample_ID` against historical records.
2. **Incremental Batch Appending**: New unique patient records are appended to the `samples` database table.
3. **Automatic Recalculation**: Upon inserting new records, the analytical engine automatically re-triggers statistical pipelines across the combined dataset, updating all overall, regional, gender, and diplotype/phenotype distributions in real time.

---

## 🔬 4. Analytical Methodology & Formulas

### 4.1 Allele & Genotype Frequencies
For a locus with Reference allele $A$ and Variant allele $a$, with valid sample size $N$:
- **Genotype Frequencies**:
  $$f(AA) = \frac{N_{AA}}{N}, \quad f(Aa) = \frac{N_{Aa}}{N}, \quad f(aa) = \frac{N_{aa}}{N}$$
- **Allele Frequencies**:
  $$p = f(A) = \frac{2 N_{AA} + N_{Aa}}{2 N}, \quad q = f(a) = \frac{2 N_{aa} + N_{Aa}}{2 N} = 1 - p$$

### 4.2 Hardy–Weinberg Equilibrium (HWE) & Chi-Square Test
- **Expected Genotype Counts**:
  $$E_{AA} = N \cdot p^2, \quad E_{Aa} = N \cdot 2pq, \quad E_{aa} = N \cdot q^2$$
- **Chi-Square Goodness-of-Fit Statistic ($\chi^2$, $df = 1$)**:
  $$\chi^2 = \frac{(N_{AA} - E_{AA})^2}{E_{AA}} + \frac{(N_{Aa} - E_{Aa})^2}{E_{Aa}} + \frac{(N_{aa} - E_{aa})^2}{E_{aa}}$$
- **P-value Calculation**:
  $$\text{P-value} = 1 - \text{CDF}_{\chi^2_1}(\chi^2)$$
- **Critical Threshold**: $\chi^2 \ge 3.841$ ($\alpha = 0.05$). Values $\ge 3.841$ indicate deviation from HWE and are highlighted with a **bright red warning badge (`⚠️ Chi² ≥ 3.841`)**.
- **Haldane Exact Test**: Automatically triggered when any expected count $E < 5$ to prevent Chi-square small-sample bias.

### 4.3 3-Tier Demographic Gazetteer Engine
Samples are assigned to 5 geographic regions using a 3-tier lookup priority:
1. **South India**: Tamil Nadu, Kerala, Karnataka, Andhra Pradesh, Telangana, Puducherry, Lakshadweep.
2. **North India**: Delhi, Punjab, Haryana, Himachal Pradesh, Jammu & Kashmir, Ladakh, Uttarakhand, Uttar Pradesh.
3. **East India**: West Bengal, Odisha, Bihar, Jharkhand, Assam, Sikkim, Arunachal Pradesh, Manipur, Meghalaya, Mizoram, Nagaland, Tripura, Andaman & Nicobar.
4. **West India**: Maharashtra, Gujarat, Goa, Rajasthan, Dadra & Nagar Haveli, Daman & Diu.
5. **Central India**: Madhya Pradesh, Chhattisgarh.

### 4.4 CPIC CYP2C19 Metabolizer Phenotype Translation Matrix

| Diplotype Call | CPIC Metabolizer Phenotype | Activity Score | Functional Description |
| :--- | :--- | :---: | :--- |
| **`*17/*17`** | **Ultrarapid Metabolizer (UM)** | $> 2.0$ | Increased enzyme activity (Hom. Gain of Function) |
| **`*1/*17`** | **Rapid Metabolizer (RM)** | $1.5$ | Increased enzyme activity (Het. Gain of Function) |
| **`*1/*1`** | **Normal Metabolizer (NM)** | $2.0$ | Normal enzyme activity (Wildtype) |
| **`*1/*2`, `*1/*3`, `*2/*17`, `*3/*17`** | **Intermediate Metabolizer (IM)** | $0.5 - 1.0$ | Reduced enzyme activity |
| **`*2/*2`, `*2/*3`, `*3/*3`** | **Poor Metabolizer (PM)** | $0.0$ | No enzyme activity (Hom. Loss of Function) |

---

## 📐 5. Project Directory Structure

```
analysis tools/
├── app.py                      # 9-Tab Streamlit Web Dashboard UI & Workflow Renderer
├── config.py                   # Marker definitions, translation rules & geographic mappings
├── core/                       # Core Analytical & Statistical Engine Package
│   ├── __init__.py
│   ├── qc.py                   # Data Quality Control, Validation & Cleaning Engine
│   ├── stats.py                # Genotype, Allele Frequency & HWE Engine
│   ├── stratification.py       # Demographic & Geographic Stratification Engine
│   └── pgx_engine.py           # Gene-agnostic rule-based diplotype/phenotype engine
├── database/                   # Data Persistence Layer
│   ├── __init__.py
│   └── db_manager.py           # SQLite database manager for batch sample growth (1k -> 10k+)
├── reports/                    # Exporter & Reporting Layer
│   ├── __init__.py
│   └── report_generator.py     # Multi-Tab Excel (.xlsx), CSV, and PDF Report Exporter
├── tests/                      # Automated Unit Test Suite (PyTest)
│   ├── __init__.py
│   ├── test_qc.py
│   ├── test_stats.py
│   ├── test_stratification.py
│   ├── test_pgx_engine.py
│   ├── test_cyp2c19.py         # Regression cases using generic rule data
│   ├── test_cyp2c9.py          # Regression cases using generic rule data
│   ├── test_db_manager.py
│   └── test_report_generator.py
├── .streamlit/                 # Streamlit UI server and theme configuration
│   └── config.toml
├── requirements.txt            # Python dependencies (Pandas, SciPy, Plotly, ReportLab, etc.)
└── README.md                   # Full platform documentation & methodology guide
```

---

## 🚀 6. How to Run & Test

### 6.1 Execute Web Application Dashboard
```bash
python3 -m streamlit run app.py
```
*Access the web dashboard in your browser at `http://localhost:8501`.*

### 6.2 Execute Automated Unit Test Suite
```bash
python3 -m pytest tests/
```
*Executes all 8 automated tests for QC, statistics, HWE Chi-Square, diplotype calling, SQLite database, and report exporters (`100% pass rate`).*

---

## 📜 7. License & Compliance
Designed for research and population pharmacogenomics analysis in compliance with CPIC guidelines.

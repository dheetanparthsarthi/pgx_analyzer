# 🔬 Complete Code Flow & Execution Trace: Upload, Mapping, Result Generation & Output

---

## 📑 Table of Contents
1. [Stage 1: File Upload & Binary Byte Ingestion](#stage-1-file-upload--binary-byte-ingestion)
2. [Stage 2: Active Dataset Decisions (Merge vs. Separate vs. Replace)](#stage-2-active-dataset-decisions-merge-vs-separate-vs-replace)
3. [Stage 3: Dynamic Column Mapping Mechanics (Auto-Detect + Manual Map)](#stage-3-dynamic-column-mapping-mechanics-auto-detect--manual-map)
4. [Stage 4: Executing the Fast Compute Pipeline (`run_fast_pipeline`)](#stage-4-executing-the-fast-compute-pipeline-run_fast_pipeline)
5. [Stage 5: Detailed Calculation Engine Flow](#stage-5-detailed-calculation-engine-flow)
   - [5.1 Quality Audit & Invalid Call Filtering (`core/qc.py`)](#51-quality-audit--invalid-call-filtering-coreqcpy)
   - [5.2 Star-Allele Calling & Phenotype Classification (`core/cyp2c19.py`)](#52-star-allele-calling--phenotype-classification-corecyp2c19py)
   - [5.3 Hardy-Weinberg Statistics & Allele Frequencies (`core/stats.py`)](#53-hardy-weinberg-statistics--allele-frequencies-corestatspy)
   - [5.4 Geographic & Demographic Stratification (`core/stratification.py`)](#54-geographic--demographic-stratification-corestratificationpy)
6. [Stage 6: Relational SQLite Storage Insertion (`database/db_manager.py`)](#stage-6-relational-sqlite-storage-insertion-databasedb_managerpy)
7. [Stage 7: Rendering the 9 Interactive Analysis Tabs (`app.py`)](#stage-7-rendering-the-9-interactive-analysis-tabs-apppy)
8. [Stage 8: Generating & Downloading PDF Audit Reports (`reports/report_generator.py`)](#stage-8-generating--downloading-pdf-audit-reports-reportsreport_generatorpy)

---

## Stage 1: File Upload & Binary Byte Ingestion

### **User Experience:**
The user navigates to Step 1 in the UI and drags & drops an Excel (`.xlsx`, `.xls`) or CSV file into the file uploader widget.

### **Responsible Code in `app.py`:**
```python
file_upload = st.file_uploader(
    "📁 Upload your patient dataset (Excel or CSV)",
    type=['csv', 'xlsx', 'xls'],
    key="file_uploader_widget"
)

if file_upload is not None:
    uploaded_file_id = f"{file_upload.name}_{file_upload.size}"
    if st.session_state.get('last_uploaded_file_id') != uploaded_file_id:
        try:
            # Step 1: Read raw binary stream from browser
            file_bytes = io.BytesIO(file_upload.getvalue())
            
            # Step 2: Parse byte stream into Pandas DataFrame
            if file_upload.name.lower().endswith('.csv'):
                df_load = pd.read_csv(file_bytes)
            else:
                df_load = pd.read_excel(file_bytes)
```

### **How it Works under the Hood:**
1. `st.file_uploader` captures the file byte payload sent over the active Streamlit WebSocket connection.
2. `file_upload.getvalue()` extracts raw binary bytes from memory.
3. `io.BytesIO(file_bytes)` converts raw bytes into a Python file-like binary stream.
4. `pd.read_csv()` or `pd.read_excel()` reads the binary stream, parses headers and data rows, and constructs a two-dimensional Pandas DataFrame (`df_load`).

---

## Stage 2: Active Dataset Decisions (Merge vs. Separate vs. Replace)

### **User Experience:**
If an active dataset is already loaded in the workspace, the system displays an alert banner with 3 decision buttons: **Merge with Active Dataset**, **Keep as Separate Dataset**, or **Replace Active Dataset**.

### **Responsible Code in `app.py`:**
```python
if has_active_data():
    st.session_state['pending_new_df'] = df_load
    st.session_state['pending_file_name'] = file_upload.name
    
if st.session_state.get('pending_new_df') is not None:
    st.info("⚠️ **An active dataset is already loaded in your workspace!** How would you like to process this new dataset file?")
    c_mrg1, c_mrg2, c_mrg3 = st.columns(3)

    # Option A: Merge Datasets
    with c_mrg1:
        if st.button("🔄 Merge with Active Dataset", key="btn_merge_dataset"):
            merged_df = pd.concat([st.session_state['uploaded_df'], st.session_state['pending_new_df']], ignore_index=True)
            clear_deletion_flag()
            st.session_state['uploaded_df'] = merged_df
            st.session_state.setdefault('datasets_dict', {})['Master Merged Dataset'] = merged_df
            st.session_state['pending_new_df'] = None
            st.session_state['workflow_step'] = '2 · Map columns'
            save_active_dataset(merged_df)
            st.rerun()

    # Option B: Keep as Separate Dataset
    with c_mrg2:
        if st.button("📁 Keep as Separate Dataset", key="btn_separate_dataset"):
            new_df = st.session_state['pending_new_df']
            fname = st.session_state.get('pending_file_name', 'New File')
            clear_deletion_flag()
            st.session_state.setdefault('datasets_dict', {})[f"Dataset 1 (Active)"] = st.session_state['uploaded_df']
            st.session_state['datasets_dict'][f"Dataset 2 ({fname})"] = new_df
            st.session_state['uploaded_df'] = new_df
            st.session_state['mapped_cols'] = {}
            st.session_state['pending_new_df'] = None
            st.session_state['workflow_step'] = '2 · Map columns'
            save_active_dataset(new_df)
            st.rerun()

    # Option C: Replace Active Dataset
    with c_mrg3:
        if st.button("🆕 Replace Active Dataset", key="btn_replace_dataset"):
            new_df = st.session_state['pending_new_df']
            fname = st.session_state.get('pending_file_name', 'New File')
            clear_deletion_flag()
            st.session_state['uploaded_df'] = new_df
            st.session_state['datasets_dict'] = {f"Active Dataset ({fname})": new_df}
            st.session_state['mapped_cols'] = {}
            st.session_state['pending_new_df'] = None
            st.session_state['workflow_step'] = '2 · Map columns'
            save_active_dataset(new_df)
            st.rerun()
```

---

## Stage 3: Dynamic Column Mapping Mechanics (Auto-Detect + Manual Map)

### **User Experience:**
In Step 2, the system automatically pre-selects matching column headers. The user reviews these pre-selections and clicks **"Confirm Column Mappings"**.

### **Responsible Code in `core/qc.py` & `app.py`:**
```python
# 1. Automatic Column Header Regex Auto-Detection in core/qc.py
def normalize_column_names(self) -> pd.DataFrame:
    df = self.raw_df.copy()
    col_mapping = {}
    
    # Apply explicit user custom mappings first
    for raw_col, target_col in self.custom_col_mapping.items():
        if raw_col in df.columns and target_col and target_col != "(None)":
            col_mapping[raw_col] = target_col
            
    # Automated Regex & Keyword Matcher
    for col in df.columns:
        if col in col_mapping: continue
        col_lower = str(col).strip().lower()
        if 'cyp2c19*2' in col_lower or 'rs4244285' in col_lower: col_mapping[col] = 'CYP2C19*2'
        elif 'cyp2c19*3' in col_lower or 'rs4986893' in col_lower: col_mapping[col] = 'CYP2C19*3'
        elif 'cyp2c19*17' in col_lower or 'rs12248560' in col_lower: col_mapping[col] = 'CYP2C19*17'
        elif 'sample' in col_lower and 'id' in col_lower: col_mapping[col] = 'Sample_ID'
        elif 'gender' in col_lower or 'sex' in col_lower: col_mapping[col] = 'Gender'
        elif 'state' in col_lower: col_mapping[col] = 'State'

    df = df.rename(columns=col_mapping)
    return df.loc[:, ~df.columns.duplicated()].copy()

# 2. Saving Mappings in app.py upon User Confirmation
if st.button("Confirm Column Mappings"):
    st.session_state['mapped_cols'] = selected_mappings
    st.session_state['workflow_step'] = '3 · Results'
    save_active_dataset(df_raw)
    st.rerun()
```

---

## Stage 4: Executing the Fast Compute Pipeline (`run_fast_pipeline`)

### **Responsible Code in `app.py`:**
```python
@st.cache_data(show_spinner=False)
def run_fast_pipeline(df: pd.DataFrame, mapped_cols_tuple: tuple):
    mapped_dict = dict(mapped_cols_tuple)
    if df is None or len(df) == 0:
        clean_df = pd.DataFrame()
        qc_report = {'total_samples': 0, ...}
        full_results = {'overall': {}, 'regional': {}, 'state_wise': {}, 'gender': {}}
        return clean_df, full_results, qc_report

    # Step A: Data Quality Control Audit
    clean_df, qc_report = DataQCEngine(df, custom_col_mapping=mapped_dict).run_qc()
    
    # Step B: Star Allele Calling & CPIC Phenotype Assignment
    pgx_df, pgx_summary = CYP2C19Translator.process_dataframe(clean_df)
    
    # Step C: Population Genetics & Demographic Stratification
    full_results = PopulationStratifier.analyze_full_dataset(pgx_df)
    
    return pgx_df, full_results, qc_report
```

---

## Stage 5: Detailed Calculation Engine Flow

### 5.1 Quality Audit & Invalid Call Filtering (`core/qc.py`)
```python
def run_qc(self) -> Tuple[pd.DataFrame, Dict[str, Any]]:
    df = self.normalize_column_names()
    
    # 1. Auto-generate Sample IDs if missing
    if 'Sample_ID' not in df.columns:
        df['Sample_ID'] = [f"SAMPLE_{i+1:05d}" for i in range(len(df))]
        
    # 2. Detect duplicate Sample IDs
    duplicated_in_file = df[df['Sample_ID'].duplicated()]['Sample_ID'].tolist()
    
    # 3. Clean Gender and State strings
    gender_map = {'male': 'Male', 'm': 'Male', 'female': 'Female', 'f': 'Female'}
    df['Gender_Clean'] = df['Gender'].astype(str).str.strip().str.lower().map(lambda x: gender_map.get(x, 'Unknown'))
    df['State_Clean'] = df['State'].astype(str).str.strip() if 'State' in df.columns else ''

    # 4. SNP Genotype Validation Loop
    invalid_genotypes_log = {}
    for snp_name, snp_info in SNP_CONFIG.items():
        if snp_name in df.columns:
            raw_vals = df[snp_name]
            
            def clean_genotype(val):
                if pd.isna(val): return np.nan
                v = str(val).strip().upper().replace('/', '')
                if v in ['AG', 'GA']: v = 'GA'
                if v in ['CT', 'TC']: v = 'CT'
                if v in snp_info['valid_genotypes']: return v
                return np.nan

            cleaned_vals = raw_vals.apply(clean_genotype)
            invalid_mask = raw_vals.notna() & cleaned_vals.isna()
            if invalid_mask.any():
                invalid_genotypes_log[snp_name] = [
                    {'sample_id': sid, 'raw_value': rval}
                    for sid, rval in zip(df.loc[invalid_mask, 'Sample_ID'], raw_vals.loc[invalid_mask])
                ]
            df[snp_name] = cleaned_vals

    return df, qc_report_dict
```

---

### 5.2 Star-Allele Calling & Phenotype Classification (`core/cyp2c19.py`)
```python
def assign_sample_diplotype(row: pd.Series) -> Tuple[str, str]:
    g_cyp2 = str(row.get('CYP2C19*2', '')).strip().upper() if pd.notna(row.get('CYP2C19*2')) else None
    g_cyp3 = str(row.get('CYP2C19*3', '')).strip().upper() if pd.notna(row.get('CYP2C19*3')) else None
    g_cyp17 = str(row.get('CYP2C19*17', '')).strip().upper() if pd.notna(row.get('CYP2C19*17')) else None

    # Count variant allele copies
    star2_count = 2 if g_cyp2 == 'AA' else (1 if g_cyp2 in ['GA', 'AG'] else 0)
    star3_count = 2 if g_cyp3 == 'AA' else (1 if g_cyp3 in ['GA', 'AG'] else 0)
    star17_count = 2 if g_cyp17 == 'TT' else (1 if g_cyp17 in ['CT', 'TC'] else 0)

    # Decision tree calling
    if star2_count == 2: diplotype = '*2/*2'
    elif star3_count == 2: diplotype = '*3/*3'
    elif star17_count == 2: diplotype = '*17/*17'
    elif star2_count == 1 and star3_count == 1: diplotype = '*2/*3'
    elif star2_count == 1 and star17_count == 1: diplotype = '*2/*17'
    elif star3_count == 1 and star17_count == 1: diplotype = '*3/*17'
    elif star2_count == 1: diplotype = '*1/*2'
    elif star3_count == 1: diplotype = '*1/*3'
    elif star17_count == 1: diplotype = '*1/*17'
    else: diplotype = '*1/*1'

    phenotype = DIPLOTYPE_PHENOTYPE_MAP.get(diplotype, 'Indeterminate')
    return diplotype, phenotype
```

---

### 5.3 Hardy-Weinberg Statistics & Allele Frequencies (`core/stats.py`)
```python
def analyze_snp(cls, series: pd.Series, snp_name: str) -> Dict[str, Any]:
    clean_series = series.dropna().astype(str).str.strip().str.upper()
    n_valid = len(clean_series)
    
    c_wt = int((clean_series == wt_gen).sum())
    c_het = int((clean_series == het_gen).sum())
    c_var = int((clean_series == var_gen).sum())
    
    # 1. Allele Frequencies
    p = ((2 * c_wt) + c_het) / (2 * n_valid)
    q = 1.0 - p
    
    # 2. Expected Genotype Counts under Equilibrium
    exp_wt = n_valid * (p ** 2)
    exp_het = n_valid * (2 * p * q)
    exp_var = n_valid * (q ** 2)
    
    # 3. Chi-Square Goodness-of-Fit Stat
    chi2_stat = (((c_wt - exp_wt) ** 2) / exp_wt) + (((c_het - exp_het) ** 2) / exp_het) + (((c_var - exp_var) ** 2) / exp_var)
    p_val = stats.chi2.sf(chi2_stat, df=1) if chi2_stat > 0 else 1.0
    exact_p_val = cls.exact_hwe_pvalue(c_wt, c_het, c_var)
    
    # 4. Threshold Boundary Check (χ² ≥ 3.4)
    final_p = exact_p_val if (exp_wt < 5 or exp_het < 5 or exp_var < 5) else p_val
    interpretation = "In HWE" if final_p >= ALPHA_SIGNIFICANCE else "Departure from HWE"
    
    return results_dict
```

---

### 5.4 Geographic & Demographic Stratification (`core/stratification.py`)
```python
def analyze_full_dataset(cls, data_processed: pd.DataFrame) -> Dict[str, Any]:
    results = {'overall': {}, 'regional': {}, 'state_wise': {}, 'gender': {}}
    
    # 1. Overall Population
    results['overall'] = analyze_subgroup(data_processed, 'Overall Population')
    
    # 2. 5 Macro-Regions
    for region in ['North India', 'South India', 'East India', 'West India', 'Central India']:
        reg_df = data_processed[data_processed['Region'] == region]
        results['regional'][region] = analyze_subgroup(reg_df, region)
        
    # 3. Individual State Profiles
    for state in data_processed['State_Clean'].unique():
        st_df = data_processed[data_processed['State_Clean'] == state]
        if len(st_df) > 0:
            results['state_wise'][str(state)] = analyze_subgroup(st_df, str(state))

    # 4. Gender Cohorts
    for gender in ['Male', 'Female']:
        gen_df = data_processed[data_processed['Gender_Clean'] == gender]
        results['gender'][gender] = analyze_subgroup(gen_df, gender)
        
    return results
```

---

## Stage 6: Relational SQLite Storage Insertion (`database/db_manager.py`)

```python
def insert_batch(self, df: pd.DataFrame) -> Tuple[int, int]:
    existing_ids = self.get_existing_sample_ids()
    records = []
    for _, row in df.iterrows():
        sid = str(row.get('Sample_ID', '')).strip()
        if not sid or sid in existing_ids: continue
        records.append((sid, row.get('Gender_Clean'), row.get('State_Clean'), ...))
        existing_ids.add(sid)

    with self._get_connection() as conn:
        cursor = conn.cursor()
        cursor.executemany("INSERT INTO samples VALUES (?, ?, ...)", records)
        conn.commit()
    return len(records), skipped_count
```

---

## Stage 7: Rendering the 9 Interactive Analysis Tabs (`app.py`)

`render_unified_results()` populates the browser screen across 9 tabs:
1. **📊 Executive Overview:** Metrics for total samples, valid counts, duplicate IDs, and invalid call warnings.
2. **🧬 Genotype & Allele:** Frequency tables and `px.bar` charts.
3. **⚖️ Hardy-Weinberg:** Observed vs. Expected counts, Chi-Square stats, and `go.Figure` line chart with $\chi^2 = 3.4$ threshold line.
4. **🌐 Population Stratification:** `px.scatter_3d` rendering 3D PCA scatter plots.
5. **💊 CYP2C19 Calling:** Metabolizer phenotype distribution matrix (PM, IM, NM, RM, UM).
6. **🛡️ Data Quality Audit:** Quality breakdown per SNP, duplicate list, and invalid call logs.
7. **🗺️ Geographic Groups:** Regional phenotype matrices for North, South, East, West, and Central India.
8. **👫 Gender-Wise:** Phenotype comparison tables for Male vs Female cohorts.
9. **📍 Dedicated State Profiles:** Selectbox `st_tab_selectbox` displaying individual state profiles.

---

## Stage 8: Generating & Downloading PDF Audit Reports (`reports/report_generator.py`)

When the user clicks **"📄 Download PDF Report"**:
```python
def generate_pdf_report(cls, full_results: dict, qc_report: dict) -> bytes:
    buffer = io.BytesIO()
    doc = SimpleDocTemplate(buffer, pagesize=letter)
    story = []
    
    story.append(Paragraph("Population Pharmacogenomics Audit Report", styles['Title']))
    story.append(Spacer(1, 12))
    # Adds executive summary tables, diplotype counts, and state matrices
    doc.build(story)
    return buffer.getvalue()
```
The resulting PDF byte stream is passed to `st.download_button()` for instant browser download.

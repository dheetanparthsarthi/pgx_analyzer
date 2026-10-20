"""
Report Exporter Engine for Multi-Tab Excel, CSV, and PDF Output.
Generates comprehensive pharmacogenomic analysis reports including:
1. Full Cleaned Dataset
2. Overall Population Summary
3. Regional Tabs (North India, South India, East India, West India)
4. State-Wise Analysis Summary
5. Gender Breakdown (Male, Female)
6. Data QC & Missing Data Audit
"""

import pandas as pd
import io
from typing import Dict, Any
import logging

from reportlab.lib.pagesizes import letter
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Table, TableStyle
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib import colors

logger = logging.getLogger(__name__)


class ReportGenerator:
    def __init__(self, full_results: Dict[str, Any] = None, qc_report: Dict[str, Any] = None):
        self.full_results = full_results or {}
        self.qc_report = qc_report or {}

    def generate_excel(self) -> bytes:
        return self.export_to_excel(self.full_results)

    def generate_pdf(self) -> bytes:
        return self.export_to_pdf(self.full_results)

    def generate_csv(self) -> bytes:
        if 'processed_dataframe' in self.full_results and isinstance(self.full_results['processed_dataframe'], pd.DataFrame) and not self.full_results['processed_dataframe'].empty:
            return self.full_results['processed_dataframe'].to_csv(index=False).encode('utf-8')
        ov_df = self._dict_to_summary_table(self.full_results.get('overall', {}))
        return ov_df.to_csv(index=False).encode('utf-8')

    @staticmethod
    def _dict_to_summary_table(stats_dict: Dict[str, Any]) -> pd.DataFrame:
        """Converts statistical dictionary into a flat Pandas DataFrame for Excel export."""
        rows = []
        
        # SNP stats
        snps = stats_dict.get('snps', {})
        for snp_name, snp_data in snps.items():
            valid_n = snp_data.get('valid_samples', 0)
            
            # Genotypes
            counts = snp_data.get('genotype_counts', {})
            pcts = snp_data.get('genotype_pcts', {})
            gen_str = ", ".join([f"{k}: {v} ({pcts.get(k, 0)}%)" for k, v in counts.items()])
            
            # Alleles
            a_counts = snp_data.get('allele_counts', {})
            a_pcts = snp_data.get('allele_pcts', {})
            allele_str = ", ".join([f"{k}: {v} ({a_pcts.get(k, 0)}%)" for k, v in a_counts.items()])
            
            # HWE
            hwe = snp_data.get('hwe', {})
            
            rows.append({
                'Category': 'SNP Analysis',
                'Parameter / Variant': snp_name,
                'Valid Samples': valid_n,
                'Genotype Distribution': gen_str,
                'Allele Distribution': allele_str,
                'Chi2 Stat': hwe.get('chi2_stat', ''),
                'P-Value': hwe.get('p_value', ''),
                'Exact P-Value': hwe.get('exact_p_value', ''),
                'HWE Interpretation': hwe.get('interpretation', '')
            })
            
        # Phenotypes
        target_gene = stats_dict.get('target_gene') or stats_dict.get('active_gene') or ''
        pgx = (
            stats_dict.get('pgx_summary')
            or stats_dict.get(f"{target_gene.lower()}_summary")
            or stats_dict.get('cyp2c19_summary')
            or stats_dict.get('cyp2c9_summary')
            or {}
        ).get('phenotypes', {})
        for pheno, pdata in pgx.items():
            rows.append({
                'Category': f"{target_gene or 'Target Gene'} Metabolizer Phenotype",
                'Parameter / Variant': pheno,
                'Valid Samples': stats_dict.get('sample_count', 0),
                'Genotype Distribution': f"Count: {pdata.get('count', 0)}",
                'Allele Distribution': f"Frequency: {pdata.get('frequency', 0)}",
                'Chi2 Stat': '-',
                'P-Value': '-',
                'Exact P-Value': '-',
                'HWE Interpretation': f"Percentage: {pdata.get('percentage', 0)}%"
            })
            
        return pd.DataFrame(rows)

    @classmethod
    def export_to_excel(cls, full_results: Dict[str, Any]) -> bytes:
        """Generates a multi-tab Excel workbook as bytes."""
        output = io.BytesIO()
        active_gene = full_results.get('active_gene') or full_results.get('target_gene') or 'CYP2C19'
        with pd.ExcelWriter(output, engine='xlsxwriter') as writer:
            # 1. Full Dataset Tab
            if 'processed_dataframe' in full_results:
                full_df = full_results['processed_dataframe']
                full_df.to_excel(writer, sheet_name='Full Cleaned Dataset', index=False)
                
            # 2. Overall Population Tab
            if 'overall' in full_results:
                ov_df = cls._dict_to_summary_table(full_results['overall'])
                ov_df.to_excel(writer, sheet_name='Overall Population', index=False)
                
            # 3. Regional Tabs (North, South, East, West)
            regional_results = full_results.get('regional', {})
            for reg_name, reg_data in regional_results.items():
                safe_reg_name = "".join([c for c in str(reg_name) if c.isalnum() or c in [' ', '_']]).strip()
                sheet_title = safe_reg_name.replace(' ', '_')[:31]  # Excel 31 char sheet limit
                if not sheet_title:
                    sheet_title = "Region"
                reg_df = cls._dict_to_summary_table(reg_data)
                reg_df.to_excel(writer, sheet_name=sheet_title, index=False)
                
            # 4. State-Wise Analysis Summary Tab
            state_results = full_results.get('state_wise', {})
            st_rows = []
            for st_name, st_data in state_results.items():
                st_snps = st_data.get('snps', {})
                lead_snp = next(iter(st_snps.keys()), None)
                lead_var_freq = 0
                if lead_snp:
                    af = st_snps[lead_snp].get('allele_freqs', {})
                    if len(af) > 1:
                        lead_var_freq = list(af.values())[1]
                    elif len(af) == 1:
                        lead_var_freq = list(af.values())[0]

                pm_data = (
                    st_data.get('pgx_summary')
                    or st_data.get(f"{active_gene.lower()}_summary")
                    or st_data.get('cyp2c19_summary')
                    or st_data.get('cyp2c9_summary')
                    or {}
                ).get('phenotypes', {})
                pm_pct = pm_data.get('Poor Metabolizer', pm_data.get('Poor Metabolizer (PM)', {})).get('percentage', 0)
                st_rows.append({
                    'State Name': st_name,
                    'Sample N': st_data.get('sample_count', 0),
                    'Missing Data Rate (%)': f"{st_data.get('missing_pct', 0)}%",
                    f"Lead Variant Freq ({lead_snp or 'N/A'})": lead_var_freq,
                    'Poor Metabolizers (%)': f"{pm_pct}%"
                })
            if st_rows:
                pd.DataFrame(st_rows).to_excel(writer, sheet_name='State-Wise Summary', index=False)

            # 5. Gender Tabs
            gender_results = full_results.get('gender', {})
            for gen_name, gen_data in gender_results.items():
                sheet_title = f"Gender_{gen_name}"[:31]
                gen_df = cls._dict_to_summary_table(gen_data)
                gen_df.to_excel(writer, sheet_name=sheet_title, index=False)

        return output.getvalue()

    @classmethod
    def export_to_pdf(cls, full_results: Dict[str, Any]) -> bytes:
        """Generates a publication-ready comprehensive PDF report."""
        buffer = io.BytesIO()
        doc = SimpleDocTemplate(buffer, pagesize=letter, rightMargin=36, leftMargin=36, topMargin=36, bottomMargin=36)
        elements = []
        
        styles = getSampleStyleSheet()
        title_style = ParagraphStyle(
            'TitleStyle',
            parent=styles['Heading1'],
            fontSize=16,
            textColor=colors.HexColor('#0F172A'),
            spaceAfter=8
        )
        heading_style = ParagraphStyle(
            'HeadingStyle',
            parent=styles['Heading2'],
            fontSize=12,
            textColor=colors.HexColor('#1E3A8A'),
            spaceBefore=12,
            spaceAfter=6
        )
        body_style = styles['BodyText']
        
        active_gene = full_results.get('active_gene') or full_results.get('target_gene') or 'CYP2C19'
        overall_stats = full_results.get('overall', {})
        total_samples = overall_stats.get('sample_count', 0)
        snps_data = overall_stats.get('snps', {})
        snp_names = list(snps_data.keys())
        snp_names_str = ", ".join(snp_names) if snp_names else "Detected Biomarkers"

        # 1. Document Header
        elements.append(Paragraph(f"Automated Population Pharmacogenomics Analysis Report: {active_gene}", title_style))
        elements.append(Paragraph(f"Development of an Automated Platform for Genotype, Allele, HWE, Diplotype, and CPIC Phenotype Analysis ({active_gene})", body_style))
        elements.append(Spacer(1, 10))
        
        # 2. Executive Metadata Summary Table
        meta_table = [
            ["Parameter", "Details"],
            ["Target Gene Analyzed", f"{active_gene}"],
            ["Total Cohort Analyzed", f"{total_samples:,} Samples"],
            ["Biomarkers / Variants", snp_names_str],
            ["Statistical Engines", "Chi-Square Goodness-of-Fit (df=1) & Haldane Exact Test"],
            ["Interpretation Framework", f"CPIC / DPWG Guidelines for {active_gene} Metabolizer Phenotypes"]
        ]
        
        t_meta = Table(meta_table, colWidths=[180, 360])
        t_meta.setStyle(TableStyle([
            ('BACKGROUND', (0, 0), (1, 0), colors.HexColor('#0F172A')),
            ('TEXTCOLOR', (0, 0), (1, 0), colors.white),
            ('FONTNAME', (0, 0), (1, 0), 'Helvetica-Bold'),
            ('BOTTOMPADDING', (0, 0), (1, 0), 5),
            ('BACKGROUND', (0, 1), (-1, -1), colors.HexColor('#F8FAFC')),
            ('GRID', (0, 0), (-1, -1), 0.5, colors.HexColor('#CBD5E1'))
        ]))
        elements.append(t_meta)
        elements.append(Spacer(1, 12))
        
        # 3. Hardy-Weinberg Equilibrium (HWE) Summary Section
        elements.append(Paragraph(f"1. Overall Population Hardy-Weinberg Equilibrium (HWE) & Allele Frequencies ({active_gene})", heading_style))
        hwe_table_data = [["Biomarker / SNP", "Valid N", "Allele Frequencies", "Chi2 (χ²)", "P-Value", "HWE Status"]]
        
        for snp_name, s_res in snps_data.items():
            hw = s_res.get('hwe', {})
            af = s_res.get('allele_freqs', {})
            freq_parts = [f"{allele}:{freq}" for allele, freq in af.items()]
            freq_str = " · ".join(freq_parts) if freq_parts else "N/A"
            hwe_table_data.append([
                snp_name,
                str(s_res.get('valid_samples', 0)),
                freq_str,
                str(hw.get('chi2_stat', 0)),
                str(hw.get('p_value', 1.0)),
                str(hw.get('interpretation', 'In HWE'))
            ])
            
        t_hwe = Table(hwe_table_data, colWidths=[100, 55, 155, 65, 65, 100])
        t_hwe.setStyle(TableStyle([
            ('BACKGROUND', (0, 0), (-1, 0), colors.HexColor('#1E3A8A')),
            ('TEXTCOLOR', (0, 0), (-1, 0), colors.white),
            ('FONTNAME', (0, 0), (-1, 0), 'Helvetica-Bold'),
            ('GRID', (0, 0), (-1, -1), 0.5, colors.HexColor('#CBD5E1'))
        ]))
        elements.append(t_hwe)
        elements.append(Spacer(1, 12))

        # 4. CPIC Diplotypes & Metabolizer Phenotypes Section
        elements.append(Paragraph(f"2. CPIC {active_gene} Diplotypes & Metabolizer Phenotypes", heading_style))
        pgx = (
            overall_stats.get('pgx_summary')
            or overall_stats.get(f"{active_gene.lower()}_summary")
            or overall_stats.get('cyp2c19_summary')
            or overall_stats.get('cyp2c9_summary')
            or {}
        )
        pheno_data = pgx.get('phenotypes', {})
        
        p_table_data = [["Metabolizer Phenotype Category", "Count (N)", "Percentage (%)"]]
        for pheno_name, pdata in pheno_data.items():
            p_table_data.append([pheno_name, str(pdata.get('count', 0)), f"{pdata.get('percentage', 0)}%"])
            
        if len(p_table_data) == 1:
            p_table_data.append(["No Phenotype Calls Recorded", "0", "0.0%"])

        t_pheno = Table(p_table_data, colWidths=[240, 130, 170])
        t_pheno.setStyle(TableStyle([
            ('BACKGROUND', (0, 0), (-1, 0), colors.HexColor('#1E40AF')),
            ('TEXTCOLOR', (0, 0), (-1, 0), colors.white),
            ('FONTNAME', (0, 0), (-1, 0), 'Helvetica-Bold'),
            ('GRID', (0, 0), (-1, -1), 0.5, colors.HexColor('#9CA3AF'))
        ]))
        elements.append(t_pheno)
        elements.append(Spacer(1, 12))

        # 5. Geographic Population Groups Section
        elements.append(Paragraph(f"3. Geographic Population Groups Matrix ({active_gene})", heading_style))
        regional_results = full_results.get('regional', {})
        
        reg_table_data = [["Geographic Region", "Sample N", "Lead Marker Freq", "Lead HWE (Chi2, p)", "Poor Metabolizer %"]]
        for r_name, r_data in regional_results.items():
            r_snps = r_data.get('snps', {})
            lead_snp = next(iter(r_snps.keys()), None)
            lead_freq_str = "N/A"
            lead_hwe_str = "N/A"
            if lead_snp:
                af = r_snps[lead_snp].get('allele_freqs', {})
                lead_freq_str = f"{lead_snp}: " + ", ".join([f"{k}={v}" for k, v in af.items()])
                hw = r_snps[lead_snp].get('hwe', {})
                lead_hwe_str = f"χ²={hw.get('chi2_stat', 0)}, p={hw.get('p_value', 1.0)}"
            
            r_phenos = (
                r_data.get('pgx_summary')
                or r_data.get(f"{active_gene.lower()}_summary")
                or r_data.get('cyp2c19_summary')
                or r_data.get('cyp2c9_summary')
                or {}
            ).get('phenotypes', {})
            pm_pct = r_phenos.get('Poor Metabolizer', r_phenos.get('Poor Metabolizer (PM)', {})).get('percentage', 0)
            
            reg_table_data.append([
                str(r_name),
                str(r_data.get('sample_count', 0)),
                lead_freq_str,
                lead_hwe_str,
                f"{pm_pct}%"
            ])
            
        if len(reg_table_data) == 1:
            reg_table_data.append(["No Regional Data Stratified", "0", "-", "-", "-"])

        t_reg = Table(reg_table_data, colWidths=[100, 60, 160, 140, 80])
        t_reg.setStyle(TableStyle([
            ('BACKGROUND', (0, 0), (-1, 0), colors.HexColor('#0F766E')),
            ('TEXTCOLOR', (0, 0), (-1, 0), colors.white),
            ('FONTNAME', (0, 0), (-1, 0), 'Helvetica-Bold'),
            ('GRID', (0, 0), (-1, -1), 0.5, colors.HexColor('#CBD5E1'))
        ]))
        elements.append(t_reg)
        elements.append(Spacer(1, 12))

        # 6. Gender Stratification Section
        elements.append(Paragraph(f"4. Gender Stratification Analysis ({active_gene})", heading_style))
        gender_results = full_results.get('gender', {})
        g_table_data = [["Cohort / Gender", "Sample N", "Normal Metabolizer %", "Intermediate %", "Poor Metabolizer %"]]
        
        for g_name in ['Male', 'Female']:
            g_data = gender_results.get(g_name, {})
            g_phenos = (
                g_data.get('pgx_summary')
                or g_data.get(f"{active_gene.lower()}_summary")
                or g_data.get('cyp2c19_summary')
                or g_data.get('cyp2c9_summary')
                or {}
            ).get('phenotypes', {})
            nm_p = g_phenos.get('Normal Metabolizer', g_phenos.get('Normal Metabolizer (NM)', {})).get('percentage', 0)
            im_p = g_phenos.get('Intermediate Metabolizer', g_phenos.get('Intermediate Metabolizer (IM)', {})).get('percentage', 0)
            pm_p = g_phenos.get('Poor Metabolizer', g_phenos.get('Poor Metabolizer (PM)', {})).get('percentage', 0)
            
            g_table_data.append([
                g_name,
                str(g_data.get('sample_count', 0)),
                f"{nm_p}%",
                f"{im_p}%",
                f"{pm_p}%"
            ])
            
        t_gen = Table(g_table_data, colWidths=[120, 70, 120, 110, 120])
        t_gen.setStyle(TableStyle([
            ('BACKGROUND', (0, 0), (-1, 0), colors.HexColor('#4338CA')),
            ('TEXTCOLOR', (0, 0), (-1, 0), colors.white),
            ('FONTNAME', (0, 0), (-1, 0), 'Helvetica-Bold'),
            ('GRID', (0, 0), (-1, -1), 0.5, colors.HexColor('#CBD5E1'))
        ]))
        elements.append(t_gen)

        doc.build(elements)
        return buffer.getvalue()

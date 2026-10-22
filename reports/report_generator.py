"""
Report Exporter Engine for Multi-Tab Excel, CSV, and Publication-Grade PDF Output.
Generates comprehensive pharmacogenomic analysis reports covering all 9 analytical views:
1. Executive Overview & Data QC Audit
2. Biomarker Genotype Distributions & Allele Frequencies
3. Hardy-Weinberg Equilibrium (HWE) Audit & Goodness-of-Fit Testing
4. Star-Allele Diplotypes & CPIC Metabolizer Phenotypes
5. Population Stratification & Subgroup Comparisons
6. Comprehensive Data Quality Control (QC) & Missing Values Audit
7. Geographic Regional Population Group Profiles (All Regions)
8. Gender-Wise Comparative Stratification (Male vs Female)
9. Dedicated State-Wise Geographic Profile Matrix (All Represented States)
"""

import pandas as pd
import io
from typing import Dict, Any, List
import logging

from reportlab.lib.pagesizes import letter
from reportlab.platypus import (
    SimpleDocTemplate, Paragraph, Spacer, Table, TableStyle, PageBreak, KeepTogether, HRFlowable
)
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib import colors
from reportlab.pdfgen import canvas

logger = logging.getLogger(__name__)


class NumberedCanvas(canvas.Canvas):
    """
    Two-pass canvas to dynamically compute total pages and add running header/footer:
    'Page X of Y' alongside document metadata.
    """
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._saved_page_states = []

    def showPage(self):
        self._saved_page_states.append(dict(self.__dict__))
        self._startPage()

    def save(self):
        num_pages = len(self._saved_page_states)
        for state in self._saved_page_states:
            self.__dict__.update(state)
            self.draw_page_decorations(num_pages)
            super().showPage()
        super().save()

    def draw_page_decorations(self, page_count: int):
        self.saveState()
        self.setFont("Helvetica", 8)
        self.setFillColor(colors.HexColor('#64748B'))

        # Running Header (pages > 1)
        if self._pageNumber > 1:
            self.drawString(36, 756, "Automated Population Pharmacogenomics Analysis Report")
            self.setStrokeColor(colors.HexColor('#E2E8F0'))
            self.setLineWidth(0.5)
            self.line(36, 750, 576, 750)

        # Running Footer (all pages)
        self.setStrokeColor(colors.HexColor('#E2E8F0'))
        self.setLineWidth(0.5)
        self.line(36, 40, 576, 40)
        self.drawString(36, 28, "Confidential · Comprehensive Multi-Tab Analytical Report · CPIC & DPWG Guidelines")
        page_text = f"Page {self._pageNumber} of {page_count}"
        self.drawRightString(576, 28, page_text)
        self.restoreState()


class ReportGenerator:
    def __init__(self, full_results: Dict[str, Any] = None, qc_report: Dict[str, Any] = None):
        self.full_results = full_results or {}
        self.qc_report = qc_report or {}
        # Ensure qc_report in full_results if available
        if self.qc_report and 'qc_report' not in self.full_results:
            self.full_results['qc_report'] = self.qc_report

    def generate_excel(self) -> bytes:
        return self.export_to_excel(self.full_results, self.qc_report)

    def generate_pdf(self) -> bytes:
        return self.export_to_pdf(self.full_results, self.qc_report)

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
            
        # Diplotypes
        target_gene = stats_dict.get('target_gene') or stats_dict.get('active_gene') or ''
        pgx_block = (
            stats_dict.get('pgx_summary')
            or stats_dict.get(f"{target_gene.lower()}_summary")
            or stats_dict.get('cyp2c19_summary')
            or stats_dict.get('cyp2c9_summary')
            or {}
        )
        diplotypes_map = pgx_block.get('diplotypes', {})
        for diplo, ddata in diplotypes_map.items():
            rows.append({
                'Category': f"{target_gene or 'Target Gene'} Diplotype Call",
                'Parameter / Variant': diplo,
                'Valid Samples': stats_dict.get('sample_count', 0),
                'Genotype Distribution': f"Count: {ddata.get('count', 0)}",
                'Allele Distribution': f"Frequency: {ddata.get('frequency', 0)}",
                'Chi2 Stat': '-',
                'P-Value': '-',
                'Exact P-Value': '-',
                'HWE Interpretation': f"Percentage: {ddata.get('percentage', 0)}%"
            })

        # Phenotypes
        pgx = pgx_block.get('phenotypes', {})
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
    def export_to_excel(cls, full_results: Dict[str, Any], qc_report: Dict[str, Any] = None) -> bytes:
        """Generates an exhaustive multi-tab Excel workbook covering all analytical options."""
        output = io.BytesIO()
        active_gene = full_results.get('active_gene') or full_results.get('target_gene') or 'CYP2C19'
        qc_data = qc_report or full_results.get('qc_report', {})

        with pd.ExcelWriter(output, engine='xlsxwriter') as writer:
            # 1. Full Dataset Tab
            if 'processed_dataframe' in full_results and isinstance(full_results['processed_dataframe'], pd.DataFrame) and not full_results['processed_dataframe'].empty:
                full_df = full_results['processed_dataframe']
                full_df.to_excel(writer, sheet_name='Full Cleaned Dataset', index=False)
                
            # 2. Overall Population Tab
            if 'overall' in full_results:
                ov_df = cls._dict_to_summary_table(full_results['overall'])
                ov_df.to_excel(writer, sheet_name='Overall Population', index=False)

            # 3. Dedicated Genotype & Allele Breakdown
            snps_data = full_results.get('overall', {}).get('snps', {})
            gt_rows = []
            for snp_name, s_res in snps_data.items():
                g_c = s_res.get('genotype_counts', {})
                g_p = s_res.get('genotype_pcts', {})
                a_c = s_res.get('allele_counts', {})
                a_p = s_res.get('allele_pcts', {})
                for gt, cnt in g_c.items():
                    gt_rows.append({
                        'Marker': snp_name,
                        'Class': 'Genotype',
                        'Call': gt,
                        'Count (N)': cnt,
                        'Percentage (%)': f"{g_p.get(gt, 0.0):.2f}%"
                    })
                for al, cnt in a_c.items():
                    gt_rows.append({
                        'Marker': snp_name,
                        'Class': 'Allele',
                        'Call': al,
                        'Count (N)': cnt,
                        'Percentage (%)': f"{a_p.get(al, 0.0):.2f}%"
                    })
            if gt_rows:
                pd.DataFrame(gt_rows).to_excel(writer, sheet_name='Genotypes & Alleles', index=False)

            # 4. Hardy-Weinberg Equilibrium (HWE) Testing
            hwe_rows = []
            overall_total = full_results.get('overall', {}).get('sample_count', 0)
            for snp_name, s_res in snps_data.items():
                hw = s_res.get('hwe', {})
                gc = s_res.get('genotype_counts', {})
                exp_c = hw.get('expected_counts', {})
                v_n = s_res.get('valid_samples', 0)
                m_n = max(0, overall_total - v_n)
                chi2_val = hw.get('chi2_stat', 0.0)
                hwe_rows.append({
                    'SNP Variant': snp_name,
                    'Cohort Total N': overall_total,
                    'Valid N': v_n,
                    'Missing N': m_n,
                    'Observed Counts': ", ".join([f"{k}:{v}" for k, v in gc.items()]),
                    'Expected Counts': ", ".join([f"{k}:{exp_c.get(k, 0)}" for k in gc.keys()]),
                    'Chi-Square (χ²)': chi2_val,
                    'df': hw.get('df', 1),
                    'P-Value': hw.get('p_value', 1.0),
                    'Exact P-Value': hw.get('exact_p_value', 1.0),
                    'HWE Status': 'DEVIATION (p < 0.05)' if chi2_val >= 3.841 else 'IN HWE'
                })
            if hwe_rows:
                pd.DataFrame(hwe_rows).to_excel(writer, sheet_name='HWE Audit', index=False)

            # 5. Dedicated Diplotype & Phenotype Analysis Tab
            ov_pgx = (
                full_results.get('overall', {}).get('pgx_summary')
                or full_results.get('overall', {}).get(f"{active_gene.lower()}_summary")
                or full_results.get('overall', {}).get('cyp2c19_summary')
                or full_results.get('overall', {}).get('cyp2c9_summary')
                or {}
            )
            diplo_rows = []
            for d_name, d_val in ov_pgx.get('diplotypes', {}).items():
                diplo_rows.append({
                    'Diplotype': d_name,
                    'Count': d_val.get('count', 0),
                    'Frequency': d_val.get('frequency', 0.0),
                    'Percentage (%)': f"{d_val.get('percentage', 0.0):.2f}%"
                })
            pheno_rows = []
            for p_name, p_val in ov_pgx.get('phenotypes', {}).items():
                pheno_rows.append({
                    'Metabolizer Phenotype': p_name,
                    'Count': p_val.get('count', 0),
                    'Frequency': p_val.get('frequency', 0.0),
                    'Percentage (%)': f"{p_val.get('percentage', 0.0):.2f}%"
                })
            if diplo_rows or pheno_rows:
                diplo_df = pd.DataFrame(diplo_rows)
                pheno_df = pd.DataFrame(pheno_rows)
                combined_pgx_df = pd.concat([diplo_df, pd.DataFrame([{} for _ in range(2)]), pheno_df], axis=1)
                combined_pgx_df.to_excel(writer, sheet_name='Diplotype & Phenotype', index=False)
                
            # 6. Regional Tabs (North, South, East, West, Central)
            regional_results = full_results.get('regional', {})
            for reg_name, reg_data in regional_results.items():
                safe_reg_name = "".join([c for c in str(reg_name) if c.isalnum() or c in [' ', '_']]).strip()
                sheet_title = safe_reg_name.replace(' ', '_')[:31]
                if not sheet_title:
                    sheet_title = "Region"
                reg_df = cls._dict_to_summary_table(reg_data)
                reg_df.to_excel(writer, sheet_name=sheet_title, index=False)
                
            # 7. State-Wise Analysis Summary Tab (All represented states)
            state_results = full_results.get('state_wise', {})
            st_rows = []
            for st_name, st_data in sorted(state_results.items(), key=lambda x: x[1].get('sample_count', 0), reverse=True):
                st_snps = st_data.get('snps', {})
                lead_snp = next(iter(st_snps.keys()), None)
                lead_var_freq = 0.0
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

            # 8. Gender Breakdown Tabs
            gender_results = full_results.get('gender', {})
            for gen_name, gen_data in gender_results.items():
                sheet_title = f"Gender_{gen_name}"[:31]
                gen_df = cls._dict_to_summary_table(gen_data)
                gen_df.to_excel(writer, sheet_name=sheet_title, index=False)

            # 9. Quality Control Audit Tab
            if qc_data:
                qc_summary_rows = [
                    {'QC Metric / Audit Parameter': 'Total Samples Processed', 'Value': qc_data.get('total_samples', 0)},
                    {'QC Metric / Audit Parameter': 'Fully Valid Genotype Samples', 'Value': qc_data.get('fully_valid_samples', 0)},
                    {'QC Metric / Audit Parameter': 'Samples With Missing Calls', 'Value': qc_data.get('samples_with_missing_count', 0)},
                    {'QC Metric / Audit Parameter': 'Overall Missing Rate (%)', 'Value': f"{qc_data.get('overall_missing_pct', 0.0):.2f}%"},
                    {'QC Metric / Audit Parameter': 'Duplicate Sample IDs', 'Value': qc_data.get('duplicate_sample_count', 0)},
                    {'QC Metric / Audit Parameter': 'Invalid Genotype Calls Detected', 'Value': qc_data.get('invalid_genotypes_count', 0)}
                ]
                pd.DataFrame(qc_summary_rows).to_excel(writer, sheet_name='QC Audit Summary', index=False)

        return output.getvalue()

    @classmethod
    def export_to_pdf(cls, full_results: Dict[str, Any], qc_report: Dict[str, Any] = None) -> bytes:
        """
        Generates a publication-ready, multi-page comprehensive PDF report.
        Exhaustively covers all 9 analytical views:
          - Page 1: Title, Executive Metadata, Executive Overview & Key Metrics (Tab 1)
          - Page 2: Biomarker Genotype Distributions & Allele Frequencies (Tab 2)
          - Page 3: Hardy-Weinberg Equilibrium (HWE) Audit & Goodness-of-Fit Testing (Tab 3)
          - Page 4: Star-Allele Diplotypes & CPIC Metabolizer Phenotypes (Tab 5)
          - Page 5: Population Stratification & Subgroup Comparisons (Tab 4)
          - Page 6: Comprehensive Data Quality Control (QC) & Missing Values Audit (Tab 6)
          - Page 7: Geographic Regional Population Group Profiles (Tab 7)
          - Page 8: Gender-Wise Comparative Stratification (Tab 8)
          - Page 9+: Dedicated State-Wise Geographic Profile Matrix (Tab 9)
        """
        buffer = io.BytesIO()
        doc = SimpleDocTemplate(
            buffer,
            pagesize=letter,
            rightMargin=36,
            leftMargin=36,
            topMargin=46,
            bottomMargin=46
        )
        elements = []
        
        # Styles
        styles = getSampleStyleSheet()
        title_style = ParagraphStyle(
            'DocTitle',
            parent=styles['Heading1'],
            fontSize=16,
            leading=20,
            textColor=colors.HexColor('#0F172A'),
            spaceAfter=4
        )
        subtitle_style = ParagraphStyle(
            'DocSubtitle',
            parent=styles['Normal'],
            fontSize=9,
            leading=12,
            textColor=colors.HexColor('#475569'),
            spaceAfter=12
        )
        tab_header_style = ParagraphStyle(
            'TabHeader',
            parent=styles['Heading2'],
            fontSize=12,
            leading=15,
            textColor=colors.HexColor('#1E3A8A'),
            spaceBefore=4,
            spaceAfter=6
        )
        section_desc_style = ParagraphStyle(
            'SectionDesc',
            parent=styles['Normal'],
            fontSize=8.5,
            leading=11.5,
            textColor=colors.HexColor('#334155'),
            spaceAfter=8
        )
        table_cell_style = ParagraphStyle(
            'TableCell',
            parent=styles['Normal'],
            fontSize=8,
            leading=10,
            textColor=colors.HexColor('#0F172A')
        )
        table_cell_bold = ParagraphStyle(
            'TableCellBold',
            parent=styles['Normal'],
            fontSize=8,
            leading=10,
            fontName='Helvetica-Bold',
            textColor=colors.HexColor('#0F172A')
        )
        table_cell_header = ParagraphStyle(
            'TableHeaderCell',
            parent=styles['Normal'],
            fontSize=8,
            leading=10,
            fontName='Helvetica-Bold',
            textColor=colors.white
        )

        active_gene = full_results.get('active_gene') or full_results.get('target_gene') or 'CYP2C19'
        overall_stats = full_results.get('overall', {})
        total_samples = overall_stats.get('sample_count', 0)
        snps_data = overall_stats.get('snps', {})
        snp_names = list(snps_data.keys())
        snp_names_str = ", ".join(snp_names) if snp_names else "Detected Biomarkers"
        qc_data = qc_report or full_results.get('qc_report', {})

        # =========================================================================
        # PAGE 1: EXECUTIVE OVERVIEW & METADATA (Tab 1)
        # =========================================================================
        elements.append(Paragraph(f"Automated Population Pharmacogenomics Analysis: {active_gene}", title_style))
        elements.append(Paragraph(f"Comprehensive Multi-Tab Clinical & Population Genetics Report · Clinical Pharmacogenetics Implementation Consortium (CPIC)", subtitle_style))
        elements.append(HRFlowable(width="100%", thickness=1, color=colors.HexColor('#0284C7'), spaceBefore=0, spaceAfter=8))
        
        # Metadata Table
        meta_table = [
            [Paragraph("Target Gene Analyzed", table_cell_bold), Paragraph(str(active_gene), table_cell_style)],
            [Paragraph("Total Cohort Analyzed", table_cell_bold), Paragraph(f"{total_samples:,} Samples", table_cell_style)],
            [Paragraph("Biomarkers / Variants", table_cell_bold), Paragraph(snp_names_str, table_cell_style)],
            [Paragraph("Statistical Engines", table_cell_bold), Paragraph("Chi-Square Goodness-of-Fit (df=1, Critical α=0.05, χ²=3.841) & Haldane Exact Test", table_cell_style)],
            [Paragraph("Clinical Guidelines", table_cell_bold), Paragraph(f"CPIC / DPWG Standardized Metabolizer Phenotype Translation ({active_gene})", table_cell_style)]
        ]
        t_meta = Table(meta_table, colWidths=[160, 380])
        t_meta.setStyle(TableStyle([
            ('BACKGROUND', (0, 0), (0, -1), colors.HexColor('#F1F5F9')),
            ('BACKGROUND', (1, 0), (1, -1), colors.HexColor('#FFFFFF')),
            ('GRID', (0, 0), (-1, -1), 0.5, colors.HexColor('#CBD5E1')),
            ('TOPPADDING', (0, 0), (-1, -1), 4),
            ('BOTTOMPADDING', (0, 0), (-1, -1), 4)
        ]))
        elements.append(t_meta)
        elements.append(Spacer(1, 12))

        # Tab 1: Executive Overview Key Metrics Table
        elements.append(Paragraph("1. Executive Overview & Cohort Quality Audit (Tab 1)", tab_header_style))
        elements.append(Paragraph("High-level audit metrics and sample validation breakdown for the analyzed cohort.", section_desc_style))

        fully_valid = qc_data.get('fully_valid_samples', total_samples)
        missing_count = qc_data.get('samples_with_missing_count', 0)
        missing_pct = qc_data.get('overall_missing_pct', 0.0)
        dup_count = qc_data.get('duplicate_sample_count', 0)
        invalid_count = qc_data.get('invalid_genotypes_count', 0)

        exec_kpi_table = [
            [Paragraph("Metric / Audit Parameter", table_cell_header), Paragraph("Value", table_cell_header), Paragraph("Percentage of Cohort", table_cell_header), Paragraph("Clinical / Audit Interpretation", table_cell_header)],
            [Paragraph("Total Samples Uploaded", table_cell_style), Paragraph(f"{total_samples:,}", table_cell_style), Paragraph("100.0%", table_cell_style), Paragraph("Total patient entries in source dataset", table_cell_style)],
            [Paragraph("Fully Valid Genotype Samples", table_cell_style), Paragraph(f"{fully_valid:,}", table_cell_style), Paragraph(f"{(fully_valid/total_samples*100):.1f}%" if total_samples>0 else "0%", table_cell_style), Paragraph("All target markers successfully genotyped", table_cell_style)],
            [Paragraph("Samples With Missing Genotypes", table_cell_style), Paragraph(f"{missing_count:,}", table_cell_style), Paragraph(f"{missing_pct:.1f}%", table_cell_style), Paragraph("Excluded from allele frequency calculations", table_cell_style)],
            [Paragraph("Duplicate Sample IDs", table_cell_style), Paragraph(f"{dup_count:,}", table_cell_style), Paragraph(f"{(dup_count/total_samples*100):.1f}%" if total_samples>0 else "0%", table_cell_style), Paragraph("Flagged in cohort for audit tracking", table_cell_style)],
            [Paragraph("Invalid / Quarantined Calls", table_cell_style), Paragraph(f"{invalid_count:,}", table_cell_style), Paragraph(f"{(invalid_count/total_samples*100):.1f}%" if total_samples>0 else "0%", table_cell_style), Paragraph("Non-canonical alleles quarantined", table_cell_style)]
        ]
        t_exec = Table(exec_kpi_table, colWidths=[150, 70, 100, 220])
        t_exec.setStyle(TableStyle([
            ('BACKGROUND', (0, 0), (-1, 0), colors.HexColor('#0F172A')),
            ('GRID', (0, 0), (-1, -1), 0.5, colors.HexColor('#CBD5E1')),
            ('ROWBACKGROUNDS', (0, 1), (-1, -1), [colors.white, colors.HexColor('#F8FAFC')]),
            ('TOPPADDING', (0, 0), (-1, -1), 4),
            ('BOTTOMPADDING', (0, 0), (-1, -1), 4)
        ]))
        elements.append(t_exec)
        elements.append(Spacer(1, 12))

        # Overall Phenotype & Regional Summary Tables Side-by-Side (Page 1)
        ov_pgx = (
            overall_stats.get('pgx_summary')
            or overall_stats.get(f"{active_gene.lower()}_summary")
            or overall_stats.get('cyp2c19_summary')
            or overall_stats.get('cyp2c9_summary')
            or {}
        )
        pheno_dict = ov_pgx.get('phenotypes', {})
        ov_p_rows = [
            [Paragraph("Metabolizer Category", table_cell_header), Paragraph("Count (N)", table_cell_header), Paragraph("Percentage (%)", table_cell_header)]
        ]
        for p_name, p_data in pheno_dict.items():
            if p_data.get('count', 0) > 0 or len(pheno_dict) <= 5:
                ov_p_rows.append([
                    Paragraph(str(p_name), table_cell_style),
                    Paragraph(str(p_data.get('count', 0)), table_cell_style),
                    Paragraph(f"{p_data.get('percentage', 0.0):.1f}%", table_cell_style)
                ])
        if len(ov_p_rows) == 1:
            ov_p_rows.append([Paragraph("No Phenotype Calls", table_cell_style), Paragraph("0", table_cell_style), Paragraph("0.0%", table_cell_style)])
        
        t_ov_p = Table(ov_p_rows, colWidths=[260, 120, 160])
        t_ov_p.setStyle(TableStyle([
            ('BACKGROUND', (0, 0), (-1, 0), colors.HexColor('#1E40AF')),
            ('GRID', (0, 0), (-1, -1), 0.5, colors.HexColor('#CBD5E1')),
            ('ROWBACKGROUNDS', (0, 1), (-1, -1), [colors.white, colors.HexColor('#F8FAFC')]),
            ('TOPPADDING', (0, 0), (-1, -1), 4),
            ('BOTTOMPADDING', (0, 0), (-1, -1), 4)
        ]))
        elements.append(Paragraph("<b>Overall Metabolizer Phenotype Distribution:</b>", section_desc_style))
        elements.append(t_ov_p)

        # =========================================================================
        # PAGE 2: GENOTYPE & ALLELE BREAKDOWN (Tab 2)
        # =========================================================================
        elements.append(PageBreak())
        elements.append(Paragraph("2. Detailed Genotype & Allele Breakdown (Tab 2)", tab_header_style))
        elements.append(Paragraph(
            "Granular observed counts, percentages, and allele frequencies for each analyzed biomarker. "
            "Calculations strictly reflect valid, non-missing patient calls.",
            section_desc_style
        ))

        genotype_table_data = [
            [
                Paragraph("Biomarker", table_cell_header),
                Paragraph("Genotype Call", table_cell_header),
                Paragraph("Zygosity Classification", table_cell_header),
                Paragraph("Count (N)", table_cell_header),
                Paragraph("Cohort Freq", table_cell_header),
                Paragraph("Cohort %", table_cell_header)
            ]
        ]
        for snp_name, s_res in snps_data.items():
            g_c = s_res.get('genotype_counts', {})
            g_p = s_res.get('genotype_pcts', {})
            g_f = s_res.get('genotype_freqs', {})
            for gt, cnt in g_c.items():
                genotype_table_data.append([
                    Paragraph(str(snp_name), table_cell_bold),
                    Paragraph(str(gt), table_cell_style),
                    Paragraph("Genotype Call", table_cell_style),
                    Paragraph(str(cnt), table_cell_style),
                    Paragraph(f"{g_f.get(gt, 0.0):.4f}", table_cell_style),
                    Paragraph(f"{g_p.get(gt, 0.0):.1f}%", table_cell_style)
                ])

        if len(genotype_table_data) == 1:
            genotype_table_data.append([Paragraph("No Data", table_cell_style), Paragraph("-", table_cell_style), Paragraph("-", table_cell_style), Paragraph("-", table_cell_style), Paragraph("-", table_cell_style), Paragraph("-", table_cell_style)])

        t_gt = Table(genotype_table_data, colWidths=[110, 80, 140, 60, 75, 75])
        t_gt.setStyle(TableStyle([
            ('BACKGROUND', (0, 0), (-1, 0), colors.HexColor('#0369A1')),
            ('GRID', (0, 0), (-1, -1), 0.5, colors.HexColor('#CBD5E1')),
            ('ROWBACKGROUNDS', (0, 1), (-1, -1), [colors.white, colors.HexColor('#F8FAFC')]),
            ('TOPPADDING', (0, 0), (-1, -1), 3),
            ('BOTTOMPADDING', (0, 0), (-1, -1), 3)
        ]))
        elements.append(t_gt)
        elements.append(Spacer(1, 10))

        # Allele Frequency Table
        elements.append(Paragraph("<b>Biomarker Allele Frequencies (2N Chromosome Allele Pool):</b>", section_desc_style))
        allele_table_data = [
            [
                Paragraph("Biomarker", table_cell_header),
                Paragraph("Allele Call", table_cell_header),
                Paragraph("Allele Type", table_cell_header),
                Paragraph("Allele Count (N)", table_cell_header),
                Paragraph("Allele Frequency (f)", table_cell_header),
                Paragraph("Frequency (%)", table_cell_header)
            ]
        ]
        for snp_name, s_res in snps_data.items():
            a_c = s_res.get('allele_counts', {})
            a_p = s_res.get('allele_pcts', {})
            a_f = s_res.get('allele_freqs', {})
            tot_a = sum(a_c.values())
            for al, cnt in a_c.items():
                freq_val = a_f.get(al, (cnt / tot_a) if tot_a > 0 else 0.0)
                allele_table_data.append([
                    Paragraph(str(snp_name), table_cell_bold),
                    Paragraph(str(al), table_cell_style),
                    Paragraph("Observed Allele", table_cell_style),
                    Paragraph(str(cnt), table_cell_style),
                    Paragraph(f"{freq_val:.4f}", table_cell_style),
                    Paragraph(f"{a_p.get(al, freq_val*100.0):.1f}%", table_cell_style)
                ])

        t_al = Table(allele_table_data, colWidths=[110, 80, 140, 60, 75, 75])
        t_al.setStyle(TableStyle([
            ('BACKGROUND', (0, 0), (-1, 0), colors.HexColor('#0F766E')),
            ('GRID', (0, 0), (-1, -1), 0.5, colors.HexColor('#CBD5E1')),
            ('ROWBACKGROUNDS', (0, 1), (-1, -1), [colors.white, colors.HexColor('#F8FAFC')]),
            ('TOPPADDING', (0, 0), (-1, -1), 3),
            ('BOTTOMPADDING', (0, 0), (-1, -1), 3)
        ]))
        elements.append(t_al)

        # =========================================================================
        # PAGE 3: HARDY-WEINBERG EQUILIBRIUM AUDIT (Tab 3)
        # =========================================================================
        elements.append(PageBreak())
        elements.append(Paragraph("3. Hardy-Weinberg Equilibrium (HWE) Audit & Goodness-of-Fit (Tab 3)", tab_header_style))
        elements.append(Paragraph(
            "Hardy-Weinberg equilibrium analysis evaluates whether genotype distributions conform to expected binomial "
            "frequencies ($p^2 + 2pq + q^2 = 1$). A calculated Chi-Square $\\chi^2 \\ge 3.841$ (or $p < 0.05, df=1$) "
            "indicates statistically significant departure from equilibrium.",
            section_desc_style
        ))

        hwe_table_data = [
            [
                Paragraph("Biomarker", table_cell_header),
                Paragraph("Valid N", table_cell_header),
                Paragraph("Observed (Hom/Het/Var)", table_cell_header),
                Paragraph("Expected (Hom/Het/Var)", table_cell_header),
                Paragraph("Chi² (χ²)", table_cell_header),
                Paragraph("df", table_cell_header),
                Paragraph("P-Value", table_cell_header),
                Paragraph("Exact P", table_cell_header),
                Paragraph("HWE Audit Status", table_cell_header)
            ]
        ]
        for snp_name, s_res in snps_data.items():
            hw = s_res.get('hwe', {})
            gc = s_res.get('genotype_counts', {})
            exp_c = hw.get('expected_counts', {})
            chi2_stat = hw.get('chi2_stat', 0.0)
            p_val = hw.get('p_value', 1.0)
            exact_p = hw.get('exact_p_value', 1.0)
            v_n = s_res.get('valid_samples', 0)

            obs_str = " / ".join([str(v) for v in gc.values()]) if gc else "-"
            exp_str = " / ".join([str(round(exp_c.get(k, 0), 1)) for k in gc.keys()]) if gc else "-"
            status_text = "⚠️ DEVIATION (χ² ≥ 3.841)" if chi2_stat >= 3.841 else "Consistent (In HWE)"

            hwe_table_data.append([
                Paragraph(str(snp_name), table_cell_bold),
                Paragraph(str(v_n), table_cell_style),
                Paragraph(obs_str, table_cell_style),
                Paragraph(exp_str, table_cell_style),
                Paragraph(f"{chi2_stat:.4f}", table_cell_style),
                Paragraph(str(hw.get('df', 1)), table_cell_style),
                Paragraph(f"{p_val:.4f}", table_cell_style),
                Paragraph(f"{exact_p:.4f}", table_cell_style),
                Paragraph(status_text, table_cell_bold if chi2_stat >= 3.841 else table_cell_style)
            ])

        t_hwe = Table(hwe_table_data, colWidths=[80, 45, 95, 95, 50, 25, 45, 45, 60])
        t_hwe.setStyle(TableStyle([
            ('BACKGROUND', (0, 0), (-1, 0), colors.HexColor('#1E3A8A')),
            ('GRID', (0, 0), (-1, -1), 0.5, colors.HexColor('#CBD5E1')),
            ('ROWBACKGROUNDS', (0, 1), (-1, -1), [colors.white, colors.HexColor('#F8FAFC')]),
            ('TOPPADDING', (0, 0), (-1, -1), 4),
            ('BOTTOMPADDING', (0, 0), (-1, -1), 4)
        ]))
        elements.append(t_hwe)
        elements.append(Spacer(1, 12))

        # Standard Chi-Square Reference Table
        elements.append(Paragraph("<b>Standard Theoretical Chi-Square Reference Distribution (df = 1):</b>", section_desc_style))
        chi2_ref_table = [
            [Paragraph("Significance Level (α)", table_cell_header), Paragraph("Critical Cutoff Value (χ², df=1)", table_cell_header), Paragraph("Interpretation & Decision Rule", table_cell_header)],
            [Paragraph("0.50 (50%)", table_cell_style), Paragraph("0.455", table_cell_style), Paragraph("Consistent with HWE (Random variation)", table_cell_style)],
            [Paragraph("0.10 (10%)", table_cell_style), Paragraph("2.706", table_cell_style), Paragraph("Non-significant variation (In HWE)", table_cell_style)],
            [Paragraph("0.05 (5%) [Standard Cutoff]", table_cell_bold), Paragraph("3.841", table_cell_bold), Paragraph("⚠️ Statistically Significant HWE Departure (p < 0.05)", table_cell_bold)],
            [Paragraph("0.01 (1%)", table_cell_style), Paragraph("6.635", table_cell_style), Paragraph("⚠️ Very Strong Departure from HWE (p < 0.01)", table_cell_style)],
            [Paragraph("0.001 (0.1%)", table_cell_style), Paragraph("10.828", table_cell_style), Paragraph("⚠️ Extremely Significant Departure (p < 0.001)", table_cell_style)]
        ]
        t_ref = Table(chi2_ref_table, colWidths=[150, 150, 240])
        t_ref.setStyle(TableStyle([
            ('BACKGROUND', (0, 0), (-1, 0), colors.HexColor('#475569')),
            ('GRID', (0, 0), (-1, -1), 0.5, colors.HexColor('#CBD5E1')),
            ('ROWBACKGROUNDS', (0, 1), (-1, -1), [colors.white, colors.HexColor('#F8FAFC')]),
            ('TOPPADDING', (0, 0), (-1, -1), 3),
            ('BOTTOMPADDING', (0, 0), (-1, -1), 3)
        ]))
        elements.append(t_ref)

        # =========================================================================
        # PAGE 4: STAR-ALLELE DIPLOTYPES & CPIC PHENOTYPES (Tab 5)
        # =========================================================================
        elements.append(PageBreak())
        elements.append(Paragraph(f"4. Star-Allele Diplotypes & CPIC Metabolizer Phenotypes: {active_gene} (Tab 5)", tab_header_style))
        elements.append(Paragraph(
            f"Pharmacogenomic translation of observed alleles into validated star (*) alleles and CPIC-defined "
            f"metabolizer phenotypes for {active_gene}.",
            section_desc_style
        ))

        # Diplotype Distribution Table
        diplo_data = ov_pgx.get('diplotypes', {})
        d_table_data = [
            [
                Paragraph("Star-Allele Diplotype Call", table_cell_header),
                Paragraph("Observed Count (N)", table_cell_header),
                Paragraph("Frequency (f)", table_cell_header),
                Paragraph("Cohort Percentage (%)", table_cell_header)
            ]
        ]
        for d_name, ddata in diplo_data.items():
            d_table_data.append([
                Paragraph(str(d_name), table_cell_bold),
                Paragraph(str(ddata.get('count', 0)), table_cell_style),
                Paragraph(f"{ddata.get('frequency', 0.0):.4f}", table_cell_style),
                Paragraph(f"{ddata.get('percentage', 0.0):.2f}%", table_cell_style)
            ])
        if len(d_table_data) == 1:
            d_table_data.append([Paragraph("No Diplotype Calls Recorded", table_cell_style), Paragraph("0", table_cell_style), Paragraph("0.0000", table_cell_style), Paragraph("0.0%", table_cell_style)])

        t_diplo = Table(d_table_data, colWidths=[180, 120, 120, 120])
        t_diplo.setStyle(TableStyle([
            ('BACKGROUND', (0, 0), (-1, 0), colors.HexColor('#0284C7')),
            ('GRID', (0, 0), (-1, -1), 0.5, colors.HexColor('#BAE6FD')),
            ('ROWBACKGROUNDS', (0, 1), (-1, -1), [colors.white, colors.HexColor('#F0F9FF')]),
            ('TOPPADDING', (0, 0), (-1, -1), 4),
            ('BOTTOMPADDING', (0, 0), (-1, -1), 4)
        ]))
        elements.append(Paragraph("<b>Diplotype Distribution Table:</b>", section_desc_style))
        elements.append(t_diplo)
        elements.append(Spacer(1, 10))

        # Metabolizer Phenotypes Table
        pheno_data = ov_pgx.get('phenotypes', {})
        p_table_data = [
            [
                Paragraph("CPIC Metabolizer Phenotype Category", table_cell_header),
                Paragraph("Observed Count (N)", table_cell_header),
                Paragraph("Frequency (f)", table_cell_header),
                Paragraph("Cohort Percentage (%)", table_cell_header)
            ]
        ]
        for pheno_name, pdata in pheno_data.items():
            p_table_data.append([
                Paragraph(str(pheno_name), table_cell_bold),
                Paragraph(str(pdata.get('count', 0)), table_cell_style),
                Paragraph(f"{pdata.get('frequency', 0.0):.4f}", table_cell_style),
                Paragraph(f"{pdata.get('percentage', 0.0):.2f}%", table_cell_style)
            ])
        if len(p_table_data) == 1:
            p_table_data.append([Paragraph("No Phenotype Calls Recorded", table_cell_style), Paragraph("0", table_cell_style), Paragraph("0.0000", table_cell_style), Paragraph("0.0%", table_cell_style)])

        t_pheno = Table(p_table_data, colWidths=[180, 120, 120, 120])
        t_pheno.setStyle(TableStyle([
            ('BACKGROUND', (0, 0), (-1, 0), colors.HexColor('#1E40AF')),
            ('GRID', (0, 0), (-1, -1), 0.5, colors.HexColor('#CBD5E1')),
            ('ROWBACKGROUNDS', (0, 1), (-1, -1), [colors.white, colors.HexColor('#F8FAFC')]),
            ('TOPPADDING', (0, 0), (-1, -1), 4),
            ('BOTTOMPADDING', (0, 0), (-1, -1), 4)
        ]))
        elements.append(Paragraph("<b>CPIC Metabolizer Classification Table:</b>", section_desc_style))
        elements.append(t_pheno)
        elements.append(Spacer(1, 10))

        # CPIC Clinical Dosing Recommendation Callout
        elements.append(Paragraph(f"<b>Standard CPIC Clinical Actionability Guidelines ({active_gene}):</b>", section_desc_style))
        if 'CYP2C9' in str(active_gene).upper():
            c_desc = (
                "<b>Warfarin:</b> Intermediate Metabolizers (*1/*2, *1/*3) require 20–50% dose reduction; "
                "Poor Metabolizers (*2/*3, *3/*3) require 50–80% reduction or DOAC alternative.<br/>"
                "<b>Phenytoin:</b> Reduce initial maintenance dose by 25% for IMs and 50% for PMs.<br/>"
                "<b>Celecoxib:</b> In PMs, reduce initial dose by 50% to prevent supratherapeutic plasma accumulation."
            )
        else:
            c_desc = (
                "<b>Clopidogrel (Plavix):</b> Intermediate (*1/*2, *1/*3) and Poor (*2/*2, *2/*3, *3/*3) metabolizers exhibit "
                "reduced active metabolite formation and higher stent thrombosis risk. Recommend alternative antiplatelet (Prasugrel/Ticagrelor).<br/>"
                "<b>Voriconazole &amp; Escitalopram:</b> In Poor Metabolizers, consider alternative antifungal (Isavuconazole) or 50% dose reduction."
            )
        rec_table = [[Paragraph(c_desc, table_cell_style)]]
        t_rec = Table(rec_table, colWidths=[540])
        t_rec.setStyle(TableStyle([
            ('BACKGROUND', (0, 0), (-1, -1), colors.HexColor('#EFF6FF')),
            ('BOX', (0, 0), (-1, -1), 1, colors.HexColor('#3B82F6')),
            ('TOPPADDING', (0, 0), (-1, -1), 6),
            ('BOTTOMPADDING', (0, 0), (-1, -1), 6),
            ('LEFTPADDING', (0, 0), (-1, -1), 8),
            ('RIGHTPADDING', (0, 0), (-1, -1), 8)
        ]))
        elements.append(t_rec)

        # =========================================================================
        # PAGE 5: POPULATION STRATIFICATION (Tab 4)
        # =========================================================================
        elements.append(PageBreak())
        elements.append(Paragraph("5. Population Stratification & Group Comparisons (Tab 4)", tab_header_style))
        elements.append(Paragraph(
            "Comparative overview of sample sizes and proportional representation across all demographic strata "
            "(Overall Cohort, Geographic Regions, and Gender Groups).",
            section_desc_style
        ))

        reg_dict = full_results.get('regional', {})
        gen_dict = full_results.get('gender', {})
        tot_n = overall_stats.get('sample_count', 0)

        group_rows = [
            [Paragraph("Population Subgroup", table_cell_header), Paragraph("Sample Size (N)", table_cell_header), Paragraph("% of Total Cohort", table_cell_header), Paragraph("Stratum Category", table_cell_header)],
            [Paragraph("Overall Analyzed Cohort", table_cell_bold), Paragraph(f"{tot_n:,}", table_cell_bold), Paragraph("100.0%", table_cell_bold), Paragraph("Full Dataset", table_cell_style)],
            [Paragraph("Region — South India", table_cell_style), Paragraph(str(reg_dict.get('South India', {}).get('sample_count', 0)), table_cell_style), Paragraph(f"{(reg_dict.get('South India', {}).get('sample_count', 0)/tot_n*100):.1f}%" if tot_n>0 else "0%", table_cell_style), Paragraph("Geographic", table_cell_style)],
            [Paragraph("Region — North India", table_cell_style), Paragraph(str(reg_dict.get('North India', {}).get('sample_count', 0)), table_cell_style), Paragraph(f"{(reg_dict.get('North India', {}).get('sample_count', 0)/tot_n*100):.1f}%" if tot_n>0 else "0%", table_cell_style), Paragraph("Geographic", table_cell_style)],
            [Paragraph("Region — East India", table_cell_style), Paragraph(str(reg_dict.get('East India', {}).get('sample_count', 0)), table_cell_style), Paragraph(f"{(reg_dict.get('East India', {}).get('sample_count', 0)/tot_n*100):.1f}%" if tot_n>0 else "0%", table_cell_style), Paragraph("Geographic", table_cell_style)],
            [Paragraph("Region — West India", table_cell_style), Paragraph(str(reg_dict.get('West India', {}).get('sample_count', 0)), table_cell_style), Paragraph(f"{(reg_dict.get('West India', {}).get('sample_count', 0)/tot_n*100):.1f}%" if tot_n>0 else "0%", table_cell_style), Paragraph("Geographic", table_cell_style)],
            [Paragraph("Region — Central India", table_cell_style), Paragraph(str(reg_dict.get('Central India', {}).get('sample_count', 0)), table_cell_style), Paragraph(f"{(reg_dict.get('Central India', {}).get('sample_count', 0)/tot_n*100):.1f}%" if tot_n>0 else "0%", table_cell_style), Paragraph("Geographic", table_cell_style)],
            [Paragraph("Gender — Female", table_cell_style), Paragraph(str(gen_dict.get('Female', {}).get('sample_count', 0)), table_cell_style), Paragraph(f"{(gen_dict.get('Female', {}).get('sample_count', 0)/tot_n*100):.1f}%" if tot_n>0 else "0%", table_cell_style), Paragraph("Sex-stratified", table_cell_style)],
            [Paragraph("Gender — Male", table_cell_style), Paragraph(str(gen_dict.get('Male', {}).get('sample_count', 0)), table_cell_style), Paragraph(f"{(gen_dict.get('Male', {}).get('sample_count', 0)/tot_n*100):.1f}%" if tot_n>0 else "0%", table_cell_style), Paragraph("Sex-stratified", table_cell_style)]
        ]

        t_grp = Table(group_rows, colWidths=[180, 110, 110, 140])
        t_grp.setStyle(TableStyle([
            ('BACKGROUND', (0, 0), (-1, 0), colors.HexColor('#0F172A')),
            ('GRID', (0, 0), (-1, -1), 0.5, colors.HexColor('#CBD5E1')),
            ('ROWBACKGROUNDS', (0, 1), (-1, -1), [colors.white, colors.HexColor('#F8FAFC')]),
            ('TOPPADDING', (0, 0), (-1, -1), 4),
            ('BOTTOMPADDING', (0, 0), (-1, -1), 4)
        ]))
        elements.append(t_grp)

        # =========================================================================
        # PAGE 6: DATA QUALITY CONTROL (QC) AUDIT (Tab 6)
        # =========================================================================
        elements.append(PageBreak())
        elements.append(Paragraph("6. Data Quality Control (QC) & Missing Values Audit (Tab 6)", tab_header_style))
        elements.append(Paragraph(
            "Detailed assessment of data integrity, genotyping call completeness, missingness rates, and duplicate screening.",
            section_desc_style
        ))

        snp_qc_stats = qc_data.get('snp_qc_stats', {})
        qc_table_data = [
            [
                Paragraph("Biomarker", table_cell_header),
                Paragraph("Valid Samples (N)", table_cell_header),
                Paragraph("Missing Calls (N)", table_cell_header),
                Paragraph("Invalid Calls (N)", table_cell_header),
                Paragraph("Missing Rate (%)", table_cell_header),
                Paragraph("Quality Audit Status", table_cell_header)
            ]
        ]
        for snp_code, sdata in snp_qc_stats.items():
            inv_cnt = sdata.get('invalid_samples', 0)
            qc_table_data.append([
                Paragraph(str(snp_code), table_cell_bold),
                Paragraph(str(sdata.get('valid_samples', 0)), table_cell_style),
                Paragraph(str(sdata.get('missing_samples', 0)), table_cell_style),
                Paragraph(str(inv_cnt), table_cell_style),
                Paragraph(f"{sdata.get('missing_pct', 0.0):.1f}%", table_cell_style),
                Paragraph("PASS" if inv_cnt == 0 else "WARNING", table_cell_bold if inv_cnt > 0 else table_cell_style)
            ])
        if len(qc_table_data) == 1:
            qc_table_data.append([Paragraph("No SNP QC Stats", table_cell_style), Paragraph("-", table_cell_style), Paragraph("-", table_cell_style), Paragraph("-", table_cell_style), Paragraph("-", table_cell_style), Paragraph("PASS", table_cell_style)])

        t_qc = Table(qc_table_data, colWidths=[130, 85, 85, 80, 80, 80])
        t_qc.setStyle(TableStyle([
            ('BACKGROUND', (0, 0), (-1, 0), colors.HexColor('#0F766E')),
            ('GRID', (0, 0), (-1, -1), 0.5, colors.HexColor('#CBD5E1')),
            ('ROWBACKGROUNDS', (0, 1), (-1, -1), [colors.white, colors.HexColor('#F8FAFC')]),
            ('TOPPADDING', (0, 0), (-1, -1), 4),
            ('BOTTOMPADDING', (0, 0), (-1, -1), 4)
        ]))
        elements.append(t_qc)
        elements.append(Spacer(1, 10))

        # Duplicate Sample IDs Audit Callout
        dup_ids = qc_data.get('duplicate_sample_ids', [])
        elements.append(Paragraph("<b>Duplicate Sample IDs Audit:</b>", section_desc_style))
        if dup_ids:
            dup_text = f"Found {len(dup_ids)} duplicate sample identifier(s): {', '.join([str(x) for x in dup_ids[:25]])}{'...' if len(dup_ids) > 25 else ''}"
        else:
            dup_text = "No duplicate sample IDs detected in uploaded cohort. All sample records possess unique identifiers."
        t_dup = Table([[Paragraph(dup_text, table_cell_style)]], colWidths=[540])
        t_dup.setStyle(TableStyle([
            ('BACKGROUND', (0, 0), (-1, -1), colors.HexColor('#F8FAFC')),
            ('BOX', (0, 0), (-1, -1), 0.5, colors.HexColor('#CBD5E1')),
            ('TOPPADDING', (0, 0), (-1, -1), 4),
            ('BOTTOMPADDING', (0, 0), (-1, -1), 4),
            ('LEFTPADDING', (0, 0), (-1, -1), 6),
            ('RIGHTPADDING', (0, 0), (-1, -1), 6)
        ]))
        elements.append(t_dup)

        # =========================================================================
        # PAGE 7: GEOGRAPHIC REGIONAL POPULATION PROFILES (Tab 7)
        # =========================================================================
        elements.append(PageBreak())
        elements.append(Paragraph(f"7. Geographic Regional Population Group Profiles (Tab 7)", tab_header_style))
        elements.append(Paragraph(
            f"Multi-region comparative distribution of sample sizes, lead variant allele frequencies, HWE test statistics, "
            f"and clinical metabolizer phenotype prevalence across Indian geographic groups.",
            section_desc_style
        ))

        regions_list = ['South India', 'North India', 'East India', 'West India', 'Central India']
        reg_table_data = [
            [
                Paragraph("Geographic Region", table_cell_header),
                Paragraph("Sample N", table_cell_header),
                Paragraph("Lead Marker Freq", table_cell_header),
                Paragraph("HWE (χ², p)", table_cell_header),
                Paragraph("Intermediate %", table_cell_header),
                Paragraph("Poor Metabolizer %", table_cell_header)
            ]
        ]
        for r_name in regions_list:
            r_data = reg_dict.get(r_name, {})
            r_snps = r_data.get('snps', {})
            lead_snp = next(iter(r_snps.keys()), None)
            lead_freq_str = "N/A"
            lead_hwe_str = "N/A"
            if lead_snp:
                af = r_snps[lead_snp].get('allele_freqs', {})
                lead_freq_str = ", ".join([f"{k}:{v}" for k, v in af.items()])
                hw = r_snps[lead_snp].get('hwe', {})
                lead_hwe_str = f"χ²={hw.get('chi2_stat', 0)}, p={hw.get('p_value', 1.0)}"
            
            r_phenos = (
                r_data.get('pgx_summary')
                or r_data.get(f"{active_gene.lower()}_summary")
                or r_data.get('cyp2c19_summary')
                or r_data.get('cyp2c9_summary')
                or {}
            ).get('phenotypes', {})
            pm_pct = r_phenos.get('Poor Metabolizer', r_phenos.get('Poor Metabolizer (PM)', {})).get('percentage', 0.0)
            im_pct = r_phenos.get('Intermediate Metabolizer', r_phenos.get('Intermediate Metabolizer (IM)', {})).get('percentage', 0.0)
            
            reg_table_data.append([
                Paragraph(str(r_name), table_cell_bold),
                Paragraph(str(r_data.get('sample_count', 0)), table_cell_style),
                Paragraph(lead_freq_str, table_cell_style),
                Paragraph(lead_hwe_str, table_cell_style),
                Paragraph(f"{im_pct:.1f}%", table_cell_style),
                Paragraph(f"{pm_pct:.1f}%", table_cell_bold if pm_pct > 0 else table_cell_style)
            ])

        t_reg = Table(reg_table_data, colWidths=[100, 60, 140, 110, 65, 65])
        t_reg.setStyle(TableStyle([
            ('BACKGROUND', (0, 0), (-1, 0), colors.HexColor('#0F766E')),
            ('GRID', (0, 0), (-1, -1), 0.5, colors.HexColor('#CBD5E1')),
            ('ROWBACKGROUNDS', (0, 1), (-1, -1), [colors.white, colors.HexColor('#F8FAFC')]),
            ('TOPPADDING', (0, 0), (-1, -1), 4),
            ('BOTTOMPADDING', (0, 0), (-1, -1), 4)
        ]))
        elements.append(t_reg)
        elements.append(Spacer(1, 10))

        # Full Regional Phenotype Distribution Matrix
        elements.append(Paragraph("<b>Regional Metabolizer Phenotype Distribution Matrix:</b>", section_desc_style))
        pheno_categories = ['Normal Metabolizer', 'Intermediate Metabolizer', 'Poor Metabolizer', 'Rapid Metabolizer', 'Ultrarapid Metabolizer']
        reg_pheno_matrix = [
            [Paragraph("Region", table_cell_header)] + [Paragraph(p.replace(' Metabolizer', ''), table_cell_header) for p in pheno_categories]
        ]
        for r_name in regions_list:
            r_data = reg_dict.get(r_name, {})
            p_dict = (
                r_data.get('pgx_summary')
                or r_data.get(f"{active_gene.lower()}_summary")
                or r_data.get('cyp2c19_summary')
                or r_data.get('cyp2c9_summary')
                or {}
            ).get('phenotypes', {})
            row = [Paragraph(str(r_name), table_cell_bold)]
            for p_cat in pheno_categories:
                p_info = p_dict.get(p_cat, p_dict.get(f"{p_cat} (PM)", p_dict.get(f"{p_cat} (IM)", p_dict.get(f"{p_cat} (NM)", {}))))
                cnt = p_info.get('count', 0)
                pct = p_info.get('percentage', 0.0)
                row.append(Paragraph(f"{cnt} ({pct:.1f}%)", table_cell_style))
            reg_pheno_matrix.append(row)

        t_reg_p = Table(reg_pheno_matrix, colWidths=[115, 85, 85, 85, 85, 85])
        t_reg_p.setStyle(TableStyle([
            ('BACKGROUND', (0, 0), (-1, 0), colors.HexColor('#1E3A8A')),
            ('GRID', (0, 0), (-1, -1), 0.5, colors.HexColor('#CBD5E1')),
            ('ROWBACKGROUNDS', (0, 1), (-1, -1), [colors.white, colors.HexColor('#F8FAFC')]),
            ('TOPPADDING', (0, 0), (-1, -1), 4),
            ('BOTTOMPADDING', (0, 0), (-1, -1), 4)
        ]))
        elements.append(t_reg_p)

        # =========================================================================
        # PAGE 8: GENDER-WISE COMPARATIVE STRATIFICATION (Tab 8)
        # =========================================================================
        elements.append(PageBreak())
        elements.append(Paragraph(f"8. Gender-Wise Comparative Stratification (Tab 8)", tab_header_style))
        elements.append(Paragraph(
            f"Sex-stratified evaluation comparing Male vs Female cohort sizes, biomarker genotype counts, and metabolizer phenotypes.",
            section_desc_style
        ))

        gender_results = full_results.get('gender', {})
        g_table_data = [
            [
                Paragraph("Cohort / Gender", table_cell_header),
                Paragraph("Sample N", table_cell_header),
                Paragraph("Normal Metabolizer %", table_cell_header),
                Paragraph("Intermediate %", table_cell_header),
                Paragraph("Poor Metabolizer %", table_cell_header),
                Paragraph("Rapid / Ultrarapid %", table_cell_header)
            ]
        ]
        for g_name in ['Male', 'Female']:
            g_data = gender_results.get(g_name, {})
            g_phenos = (
                g_data.get('pgx_summary')
                or g_data.get(f"{active_gene.lower()}_summary")
                or g_data.get('cyp2c19_summary')
                or g_data.get('cyp2c9_summary')
                or {}
            ).get('phenotypes', {})
            nm_p = g_phenos.get('Normal Metabolizer', g_phenos.get('Normal Metabolizer (NM)', {})).get('percentage', 0.0)
            im_p = g_phenos.get('Intermediate Metabolizer', g_phenos.get('Intermediate Metabolizer (IM)', {})).get('percentage', 0.0)
            pm_p = g_phenos.get('Poor Metabolizer', g_phenos.get('Poor Metabolizer (PM)', {})).get('percentage', 0.0)
            rm_p = g_phenos.get('Rapid Metabolizer', {}).get('percentage', 0.0) + g_phenos.get('Ultrarapid Metabolizer', {}).get('percentage', 0.0)
            
            g_table_data.append([
                Paragraph(str(g_name), table_cell_bold),
                Paragraph(str(g_data.get('sample_count', 0)), table_cell_style),
                Paragraph(f"{nm_p:.1f}%", table_cell_style),
                Paragraph(f"{im_p:.1f}%", table_cell_style),
                Paragraph(f"{pm_p:.1f}%", table_cell_style),
                Paragraph(f"{rm_p:.1f}%", table_cell_style)
            ])

        t_gen = Table(g_table_data, colWidths=[100, 70, 95, 95, 90, 90])
        t_gen.setStyle(TableStyle([
            ('BACKGROUND', (0, 0), (-1, 0), colors.HexColor('#4338CA')),
            ('GRID', (0, 0), (-1, -1), 0.5, colors.HexColor('#CBD5E1')),
            ('ROWBACKGROUNDS', (0, 1), (-1, -1), [colors.white, colors.HexColor('#F8FAFC')]),
            ('TOPPADDING', (0, 0), (-1, -1), 4),
            ('BOTTOMPADDING', (0, 0), (-1, -1), 4)
        ]))
        elements.append(t_gen)
        elements.append(Spacer(1, 10))

        # Gender Genotype Comparison Table for All Biomarkers
        elements.append(Paragraph("<b>Biomarker Genotype Distributions (Male vs Female):</b>", section_desc_style))
        gen_gt_table = [
            [
                Paragraph("Biomarker", table_cell_header),
                Paragraph("Genotype", table_cell_header),
                Paragraph("Male Count (N)", table_cell_header),
                Paragraph("Male %", table_cell_header),
                Paragraph("Female Count (N)", table_cell_header),
                Paragraph("Female %", table_cell_header)
            ]
        ]
        male_data = gender_results.get('Male', {})
        female_data = gender_results.get('Female', {})
        for snp_code in snps_data.keys():
            m_gts = male_data.get('snps', {}).get(snp_code, {}).get('genotype_counts', {})
            m_pcts = male_data.get('snps', {}).get(snp_code, {}).get('genotype_pcts', {})
            f_gts = female_data.get('snps', {}).get(snp_code, {}).get('genotype_counts', {})
            f_pcts = female_data.get('snps', {}).get(snp_code, {}).get('genotype_pcts', {})
            all_gts = sorted(list(set(list(m_gts.keys()) + list(f_gts.keys()))))
            for gt in all_gts:
                gen_gt_table.append([
                    Paragraph(str(snp_code), table_cell_bold),
                    Paragraph(str(gt), table_cell_style),
                    Paragraph(str(m_gts.get(gt, 0)), table_cell_style),
                    Paragraph(f"{m_pcts.get(gt, 0.0):.1f}%", table_cell_style),
                    Paragraph(str(f_gts.get(gt, 0)), table_cell_style),
                    Paragraph(f"{f_pcts.get(gt, 0.0):.1f}%", table_cell_style)
                ])

        if len(gen_gt_table) == 1:
            gen_gt_table.append([Paragraph("No Data", table_cell_style), Paragraph("-", table_cell_style), Paragraph("-", table_cell_style), Paragraph("-", table_cell_style), Paragraph("-", table_cell_style), Paragraph("-", table_cell_style)])

        t_gen_gt = Table(gen_gt_table, colWidths=[120, 80, 85, 85, 85, 85])
        t_gen_gt.setStyle(TableStyle([
            ('BACKGROUND', (0, 0), (-1, 0), colors.HexColor('#6D28D9')),
            ('GRID', (0, 0), (-1, -1), 0.5, colors.HexColor('#CBD5E1')),
            ('ROWBACKGROUNDS', (0, 1), (-1, -1), [colors.white, colors.HexColor('#F8FAFC')]),
            ('TOPPADDING', (0, 0), (-1, -1), 3),
            ('BOTTOMPADDING', (0, 0), (-1, -1), 3)
        ]))
        elements.append(t_gen_gt)

        # =========================================================================
        # PAGE 9+: DEDICATED STATE-WISE GEOGRAPHIC PROFILES (Tab 9)
        # =========================================================================
        state_results = full_results.get('state_wise', {})
        if state_results:
            elements.append(PageBreak())
            elements.append(Paragraph(f"9. Dedicated Individual State Profiles Matrix (Tab 9)", tab_header_style))
            elements.append(Paragraph(
                f"Exhaustive geographic summary across all individual states present in the cohort. "
                f"Displays sample counts (N), missing data rate, lead variant allele frequency, and poor metabolizer percentage.",
                section_desc_style
            ))

            st_table_data = [
                [
                    Paragraph("State Name", table_cell_header),
                    Paragraph("Sample N", table_cell_header),
                    Paragraph("Missing %", table_cell_header),
                    Paragraph("Lead Variant Freq", table_cell_header),
                    Paragraph("Poor Metabolizer %", table_cell_header)
                ]
            ]
            for st_name, st_data in sorted(state_results.items(), key=lambda x: x[1].get('sample_count', 0), reverse=True):
                st_snps = st_data.get('snps', {})
                lead_snp = next(iter(st_snps.keys()), None)
                lead_var_freq = 0.0
                if lead_snp:
                    af = st_snps[lead_snp].get('allele_freqs', {})
                    if len(af) > 1:
                        lead_var_freq = list(af.values())[1]
                    elif len(af) == 1:
                        lead_var_freq = list(af.values())[0]

                st_pm = (
                    st_data.get('pgx_summary')
                    or st_data.get(f"{active_gene.lower()}_summary")
                    or st_data.get('cyp2c19_summary')
                    or st_data.get('cyp2c9_summary')
                    or {}
                ).get('phenotypes', {})
                st_pm_pct = st_pm.get('Poor Metabolizer', st_pm.get('Poor Metabolizer (PM)', {})).get('percentage', 0.0)

                st_table_data.append([
                    Paragraph(str(st_name), table_cell_bold),
                    Paragraph(str(st_data.get('sample_count', 0)), table_cell_style),
                    Paragraph(f"{st_data.get('missing_pct', 0.0):.1f}%", table_cell_style),
                    Paragraph(f"{lead_var_freq:.4f}", table_cell_style),
                    Paragraph(f"{st_pm_pct:.1f}%", table_cell_bold if st_pm_pct > 0 else table_cell_style)
                ])

            t_state = Table(st_table_data, colWidths=[160, 80, 90, 110, 100])
            t_state.setStyle(TableStyle([
                ('BACKGROUND', (0, 0), (-1, 0), colors.HexColor('#B45309')),
                ('GRID', (0, 0), (-1, -1), 0.5, colors.HexColor('#FDE68A')),
                ('ROWBACKGROUNDS', (0, 1), (-1, -1), [colors.white, colors.HexColor('#FFFBEB')]),
                ('TOPPADDING', (0, 0), (-1, -1), 3),
                ('BOTTOMPADDING', (0, 0), (-1, -1), 3)
            ]))
            elements.append(t_state)

        # Build Document using NumberedCanvas for dynamic page numbering
        doc.build(elements, canvasmaker=NumberedCanvas)
        return buffer.getvalue()

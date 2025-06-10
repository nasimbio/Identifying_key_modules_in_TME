# identifying key modules that responsible for intrinsic and extrinsic mechanisms of tumor development in tumor micro-environment (TME), from 21 mouse scRNA libraries.

Project Overview
This project analyzes 21 single-nuclei RNA-seq samples across three tumor cell lines, each treated with either SOS or Veh treatments. 
The goal is to understand treatment effects on cell composition, tumor pathway activity, and tumor-specific gene expression patterns.
---

## 📊 Analysis Tasks

### **Task 1: Initial Processing and Broad Annotation**
- QC, filtering, normalization, clustering
- Broad cell type annotation
- Composition comparison across treatments and cell lines

### **Task 2: Tumor vs. Normal Cell Identification**
- YFP expression mapped via Kallisto
- InferCNV CNV profiling
- Tumor cell flagging for downstream analysis

### **Task 3: DEG and GO Enrichment in Epithelial Cells**
- Single-cell and pseudobulk DGE
- DEG overlap across 3 cell lines
- GO enrichment with clusterProfiler

### **Task 4: Pathway Scoring**
- Scoring of hallmark gene sets per cell
- Pathway trends compared across tumors and treatments

### **Task 5: COX Pathway Gene Expression in Epithelial cells**
- Identify Ptges+, Ptgs1+, Ptgs2+ cells
- Analyze composition changes by treatment and cell line
- Co-expression statistics and visualizations

---

## 💡 Notes

- Raw data is not publicly available due to client ownership and confidentiality.
- Some example outputs plots are organized by task in the `output/` folder.
- This project is designed for both reproducibility and clarity.

---

## 📬 Contact

*Author:* Nasim Rahmatpour 
*Email:* nasimrahmatpour1@gmail.com 
*GitHub:* (https://github.com/nasimbio)


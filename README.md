# Identifying key modules that are responsible for intrinsic and extrinsic mechanisms of tumor development in tumor micro-environment (TME), from 21 mouse scRNA libraries.

Project Overview
This project analyzes 21 single-nuclei RNA-seq samples across three tumor cell lines, each treated with either SOS or Veh treatments. 
The goal is to understand treatment effects on cell composition, tumor pathway activity, and tumor-specific gene expression patterns.
---

## 📊 Analysis Tasks

### **Task 1: Initial Processing and Broad Annotation**
- QC, filtering, normalization, clustering
- Broad cell type annotation
- Cell type composition comparison across treatments and cell lines

### **Task 2: Tumor vs. Normal Cell Identification**
- Quantifying the expression of YFP gene (The Yellow Florescence gene was introduced into the tumor cells, and its expression is used as a marker to identify and track tumor cells)
- Copy Number Variation Analysis (InferCNV)


### **Task 3: DEG and GO Enrichment in Epithelial Cells (how the genes are differentially expressed in response to the treatment at different cell line)**
- DGE & GO Enrichment Analysis (at single cell level)
- DGE & GO Enrichment Analysis (pseudo bulk analysis)
- Finding DEGs overlap across 3 cell lines

### **Task 4: Pathway Scoring and Visualization on Epithelial/tumor cells**
- Scoring of hallmark gene sets per cell
- Comparing the Pathway trends in tumor cells, across treatments in different cell lines

### **Task 5: COX Pathway Gene Expression in Epithelial cells**
- Identifying the cells express the COX pathway genes (Ptges+, Ptgs1+, Ptgs2+ cells)
- Analyzing the co-expression of the genes, providing statistics and visualizations

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


# design 
## williams 2023


Data from or derived Williams et al., 2023.

### Most tables of interest

```mermaid 
erDiagram

rLogCountsGFPplus {
int WBINT PK
char WBID 
float embryo_rep1 
float embryo_rep2 
float embryo_rep3 
float L1_rep1 
float L1_rep3 
float L3_rep1 
float L3_rep2 
float L3_rep3 
float embryo 
float L1 
float L3 
timestamp update_time 
}

rlogCountsLong {
  char WBID 
  varchar GENE_NAME
  int WBINT
  varchar value
  varchar stage
  int replicate 
  varchar sampleSource
  varchar geneName
  timestamp update_time
    }

allCounts {
  int allCountsID PK
  int WBINT
  char WBID
  float value
  varchar stage
  varchar sampleSource
  int replicate
  varchar transformation
  timestamp update_time 
}

log2FoldChangeWide {
int log2FoldChange_pk PK
char WBID 
varchar geneName 
float baseMean 
float log2FoldChange 
float lfcSE 
float stat_enriched 
float stat_depleted 
float stat_equal 
float pvalue_enriched 
float pvalue_depleted 
float pvalue_equal 
float padj_enriched 
float padj_depleted 
float padj_equal 
varchar outcome_01 
varchar outcome_05 
varchar stage 
timestamp update_time
}


intestineGeneCategories { 
int idintestineGeneCategories PK
char WBID 
varchar altHypDESeq 
varchar intestineExpression 
varchar stage 
float significanceCutoff 
int WBINT 
timestamp update_time
}

```

### Data from Dineen et al., that Rob reanalyzed

```mermaid

erDiagram

dineenSetsAnalyzed {
  int(11) WBINT  PK
  char(14) WBID
  float baseMean 
  float log2FoldChange 
  float lfcSE 
  float stat 
  float pvalue 
  float padj 
  varchar(45) wormbase_gseq 
  varchar(45) wikigene_name 
  varchar(45) status 
  varchar(45) description
  timestamp update_time 
}

dineenSourceData {
int WBINT PK
char WBID 
int wt_sorted_1 
int wt_sorted_2 
int wt_sorted_3 
int wt_sorted_4 
int elt7D_sorted_1 
int elt7D_sorted_2 
int elt7D_sorted_3 
int elt2D_sorted_1 
int elt2D_sorted_2 
int elt2D_sorted_3 
int elt2D_sorted_4 
int elt2Delt7D_sorted_1 
int elt2Delt7D_sorted_2 
int elt2Delt7D_sorted_3 
timestamp update_time
}

```

### Rob's raw data

All one table (rawCounts). The displayed '_stage' labels are for the graph display only. Need to check up on embryo rep1 whole and rep3 cells.

```mermaid
erDiagram

rawCounts {
int countsId PK
int WBINT
varchar WBID
float countsPerMillion
timestamp update_time
}

rawCounts_embryo {
int embryo_cells_rep1
int embryo_GFPplus_rep1
int embryo_GFPminus_rep1
int embryo_whole_rep2
int embryo_cells_rep2
int embryo_GFPplus_rep2
int embryo_GFPminus_rep2
int embryo_whole_rep3
int embryo_GFPplus_rep3
int embryo_GFPminus_rep3
}

rawCounts_L1 {
int L1_whole_rep1
int L1_cells_rep1
int L1_GFPplus_rep1
int L1_GFPminus_rep1
int L1_whole_rep2
int L1_cells_rep2
int L1_GFPplus_rep2
int L1_GFPminus_rep2
int L1_whole_rep3
int L1_cells_rep3
int L1_GFPplus_rep3
int L1_GFPminus_rep3
}

rawCounts_L3 {
int L3_whole_rep1
int L3_cells_rep1
int L3_GFPplus_rep1
int L3_GFPminus_rep1
int L3_whole_rep2
int L3_cells_rep2
int L3_GFPminus_rep2
int L3_GFPplus_rep2
int L3_whole_rep3
int L3_cells_rep3
int L3_GFPplus_rep3
int L3_GFPminus_rep3
}
```



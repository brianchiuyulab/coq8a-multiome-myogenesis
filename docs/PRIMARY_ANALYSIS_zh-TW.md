# 從頭、固定條件的主分析：先看全貌，再談位點

這份文件把「不為了得到 1.28 而調條件」作為主分析的原則。這是已看過資料後重新固定的分析路徑，**不是事前註冊**。所稱「無偏差」是指在事先指定的肌肉分化範圍內，讓所有合格基因與 peak 進入同一套檢驗和多重校正；它不是全基因組所有生物路徑的無假說探索。

## 在看結果前固定的主流程

1. **問題範圍**：MSigDB Hallmark Myogenesis ∪ Reactome Myogenesis，在此資料可測量的 221 個基因。MYOD1、MYOG、MYF5 都在母集合裡，沒有先選 MYOD1 或 MRF motif。
2. **同核資料與 QC**：四個 GSE208248 10x Multiome library。RNA UMI≥500、開放 ATAC peak≥500、各 library 去除兩種深度最高 5%、FRiP≥0.25、nucleosome signal<4、blacklist 比例<0.05、TSS enrichment≥3。QC 在看候選位點效果前完成。
3. **主要 COQ8A 比較**：每個 library 內，COQ8A≥2 UMI 與恰好 1 UMI 的核，以 RNA UMI 與 ATAC 開放 peak 數做 1:1 配對；共 958 對。≥3 UMI 與 TSS≥2 只列敏感度，不用來決定主結果。
4. **候選區域**：221 基因任一蛋白編碼 transcript TSS 的 ±100 kb、至少三個 library 有對應的 peak；去 blacklist、限制寬度與染色體、以座標做 library 間對應。共 7,699 個 peak–gene 對。納入時不看 COQ8A 差異或 motif 分數。
5. **先報全域結果**：兩個原先選定的 RNA programme、附近 ATAC programme、所有 221 個基因的 ATAC 效果，以及所有可檢驗的候選 peak。各完整檢驗集合各自做 BH 多重校正。只有看完這一步，才討論排名和任何特定位點。

程式入口是 [`run_all.py`](../run_all.py)，主結果逐項列在 [`primary_decision_summary.tsv`](../results/tables/primary_decision_summary.tsv)，全體 peak 效果與 q 值在 [`candidate_peak_effects_pooled.tsv.gz`](../results/tables/candidate_peak_effects_pooled.tsv.gz)。[Figure 1](../figures/main/Figure_1_global_discovery.pdf) 顯示設計、RNA 與 ATAC programme 和全部 221 基因；[全候選 peak 補圖](../figures/supplement/Supplementary_Figure_Global_ATAC_Scan.pdf) 顯示完整 ATAC 檢驗分布。

## 主分析實際得到的答案

| 項目 | 主結果 |
|---|---|
| RNA Hallmark Myogenesis | 4/4 library 同方向；pair-level p=0.02684，四 programme q=0.03579 |
| ATAC Hallmark Myogenesis 附近區域 | 2/4 同方向；p=0.05187，四 programme q=0.20747 |
| 221 個基因附近 ATAC 區域 | 0 個 q<0.05；最小 q≈0.302 |
| 可檢驗的候選 ATAC peaks | 6,871 個；0 個 q<0.05；最小 q≈0.119 |
| MYOD1 19 個跨四 library 共同 peak | 1.0688 倍，+0.824 個百分點；p=0.03952，但 221 基因 q=0.63537 |

**因此，這條固定主流程會得到「COQ8A 與肌肉分化 RNA programme 有關聯」，但不會得到通過完整多重校正的染色質 target。** 它不會自行導向「MYOD1 六區域 1.28 倍」作為主要發現。

## 1.28 放在哪裡

在完整主結果之後，MYOD1 可作為探索性位點追查。用 TSS≥2 核建立 peak–RNA links，選出六個低 MRF motif 分數的 MYOD1 linked peaks，再於 TSS≥3 的 COQ8A≥3 對 1 UMI 的 201 對核檢驗，**1.280 倍、+2.902 個百分點、p=0.02558、4/4 同方向**可精確重現。若建立 links 時也要求 TSS≥3，其中一個 peak 在第四個 library 只有 8 個開放核，低於 10 個的 link 檢驗門檻，於是變成五個 peak、1.241 倍。差異是**位點集合的建立條件**，不是核配對或 fold 計法被偷偷更換。[完整對照](../results/tables/myod1_link_gate_reconciliation.tsv)與[20 條敏感度結果](../results/tables/myod1_locus_summary.tsv)都保留。

1.28 可以是後續實驗的探索線索；它沒有改變主分析「沒有 FDR 通過的 ATAC target」的答案。這樣寫可以同時忠實呈現你關心的效果與完整搜尋的結果。

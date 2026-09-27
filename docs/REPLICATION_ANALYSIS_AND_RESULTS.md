# CAV3 prioritization and GSE240061 sensitivity analysis

Completed analysis, 2026-09-27. This report supersedes the *analysis-status*
statements in the earlier feasibility report; it does not change its input
inventory or the previously reported GSE208248 calculations.

The subsequent [local RNA-target follow-up](DIFFERENTIATION_SCOPE_AND_CHR16_RNA.md)
tested the chr16 interval against available neighboring RNAs. Its UBE2I link
was near zero and UBE2I RNA did not increase in the matched COQ8A-high group;
the interval is consequently deprioritized as a differentiation-mechanism lead.

## 結論先讀

1. **可以用外部功能資料合理地把 CAV3 放進候選名單。**最直接的路線是
   myogenesis → fusion 功能 → promoter 位置 → 外部肌肉調控註記。
   在本次探索的較窄範圍，CAV3 與 KCNH1 一起留下，並非只測 CAV3。
2. **舊資料有一組 CAV3 FC=1.824、p=0.0385、q=0.0770 的設定。**
   但同一候選集合改 TSS≥3 後 q=0.1065；擴大 promoter 窗口後 q 也上升。
   適合作為候選優先順序，不宣稱這條路是最初盲選得到的唯一答案。
3. **新資料確實能找到 CAV3 的有利結果。**最低 p 的設定得到 FC=3.194、
   p=0.000860、fusion-promoter q=0.01376；有另一設定 FC=3.628、q=0.0221。
   這些大倍率都出現在未配對、加入 RNA 深度共變項的模型。
4. **配對後最值得保留的 CAV3 線索約 1.32 倍。**六位 donor 的開放比例
   都增加，精確 donor sign-flip p=0.03125；計數模型 FC=1.329、p=0.0849、
   q=0.3678。換寬鬆配對又降到約 1.05–1.06 倍。
5. **目前沒有得到跨設定穩定、且 RNA 同步增加的 CAV3 機制鏈。**可列為
   探索候選，尚不值得只憑「3 倍」把整個驗證計畫押在 CAV3。
6. 另有 **chr16:1311478–1312392**，在配對後 FC≈3、完整候選範圍
   q≈0.02。它因落在 CACNA1H 的100-kb窗口進入候選，但補查全基因註記後，
   實際落在 **UBE2I 基因內並覆蓋一個 UBE2I transcript TSS**。不能將其
   當作已配對 CACNA1H 的調控位點；target RNA 尚未確立。

這次共完成舊資料 **288 組候選範圍／設定比較**、新資料 **448 組設定**。
新資料 393 組可估 ATAC/RNA 計數模型；55 組因 donor、設計自由度或計數不足
未能估計。CAV3 本身在 358 組通過計數可評估條件。未能估計不記成「沒有差」。

## 1. 外部依據如何導向 CAV3

基因全集沿用已凍結的 Hallmark/Reactome myogenesis 221 genes；
GSE208248 四個 library 共同候選 peaks 為 5,097。
外部 GO fusion 註記與 221 genes 交集為 12 genes：

ADAM12, CACNA1H, CAV3, CDON, ITGB1, KCNH1, MAPK14, MYH9, MYOD1,
MYOG, NEO1, NOS1。

這條功能路線不要求先屬早／中／晚期。以既有轉錄起始位置註記取
peak midpoint 到相應 gene TSS ≤500 bp，得到 **21 個 distinct peaks**。
用已凍結的 HSMM strong-enhancer annotation 加上八種非肌肉參考細胞皆無
strong-enhancer annotation 的條件，留下 **2 peaks**：

| Gene | hg38 interval | COQ8A-high open / low open | FC | p | q |
|---|---|---:|---:|---:|---:|
| CAV3 | chr3:8733438–8733968 | 31 / 17 | 1.8235 | 0.038477 | 0.076955 |
| KCNH1 | chr1:211133744–211134598 | 68 / 63 | 1.0794 | 0.660884 | 0.660884 |

以上為 **TSS≥2、COQ8A≥3 vs exactly 1 UMI、212 對核**；各 library 內沿用
雙深度 matching。生物來源仍為 **2 個 donor-derived lines、4 libraries**。
表中 p 是 matched-nucleus discordance 的雙側 binomial test；q 在這兩個
distinct peaks 內做 BH，不是 donor-level p，也不是全 5,097 peaks 的 q。

這是**對已知候選進行外部依據與參數的回溯探索**。外部 memberships 不使用
COQ8A 結果，但選擇哪種 membership/window 作重點展示已看過結果。
所有 288 組結果保留，不把這次路線改寫成原始事前計畫。

### 同一邏輯的鄰近設定

| Change, keeping COQ8A≥3 vs 1 | CAV3 FC | CAV3 p | q in stated scope |
|---|---:|---:|---:|
| Fusion, ±500-bp midpoint, muscle-selective; TSS≥2 | 1.824 | 0.0385 | 0.0770; 2 peaks |
| Same scope; TSS≥3 | 1.765 | 0.0533 | 0.1065; 2 peaks |
| Fusion, ±2-kb midpoint, muscle-selective; TSS≥2 | 1.824 | 0.0385 | 0.1154; 3 peaks |
| All fusion promoters within 500 bp; TSS≥2 | 1.824 | 0.0385 | 0.2693; 21 peaks |

500 bp 的意義是聚焦 TSS 近端、減少 gene assignment 距離的不確定；
肌肉選擇性的意義是聚焦外部參考中較偏肌肉的調控區域。
兩項都有生物理由，但都不是所有 fusion 調控區域的必要條件，也不能用
「有理由」代替參數敏感度結果。圖 1 顯示全部 21 peaks 的外部註記，
讓最後留下兩個位點的原因可以直接核查。

CAV3 focal midpoint 距最近註記 TSS 99 bp，區間包含該 TSS。這個 promoter
位置與另一套 chromatin-state 註記為 strong-enhancer 可以並存；兩者不是
實驗證明此 peak 必然調控 CAV3 的替代品。MYOD1 occupancy 顯示於圖 1，
但**不是**兩位點路線新增的必要篩選條件。

## 2. 新資料的樣本、前處理與統計

### Study and input

[GSE240061](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE240061)
對應 [原研究](https://pmc.ncbi.nlm.nih.gov/articles/PMC12212352/)：
6 位受試者，4 位 exercise（E/G/I/J）、2 位 rest（L/N），每人前後採樣，
12 biopsies。exercise 為 40-min cycling，post biopsy 為其後 3.5 h。
原文 Discussion 明述六人皆為 young healthy adults；沒有 young/old 分組。
每顆核同時有 RNA/ATAC；作者保留 37,154 nuclei、144,663 ATAC peaks。
保留原作者 14 類 `refined_annotations_wknn_0.8` 註記。

分析分開 Satellite Cells、Fast、Slow、Intermediate，並分開 Pre/Post。
**沒有把前後採樣當 12 個獨立 donor**；每個時間點內以 donor block 比高低。
exercise/rest 分組另作選定候選的敏感度；沒有估計 exercise × time 的因果效應。

作者 RDS 沒有可用的每核 TSS/FRiP/nucleosome metrics。因此本次使用作者
保留的核，另外測試每個 sample×celltype 的 RNA depth 與 detected-ATAC
depth 各自 ≤95th percentile 的 cap。這是 depth sensitivity，**不是新算 TSS QC**。

### Frozen coordinate mapping

舊 peaks 與新作者 peaks 要求雙向各 ≥50% overlap。不同 scopes 映射為：

| Old scope | Distinct new author peaks |
|---|---:|
| 5,097 shared myogenesis candidates | 3,071 |
| 1,777 HSMM regulatory peaks | 1,026 |
| 401 muscle-selective peaks | 191 |
| 31 fusion promoters within 2 kb | 23 |

新 CAV3 peak 為 **chr3:8733222–8734203**，包含舊 530-bp peak，
但新 peak 寬 981 bp；兩次計數窗口並非完全一致。未映射到作者 peak
代表此候選在指定對應规则下不可評估，不代表染色質關閉。

### Complete grid

4 cell types × 2 time points × 2 depth-QC settings × 7 COQ8A definitions ×
4 matching settings = **448 configurations**。

- Raw UMI definitions: ≥2 vs 1, ≥3 vs 1, ≥4 vs 1, ≥3 vs 1–2。
- Positive-only quantiles: 在每個 sample×celltype×QC 內、COQ8A>0 的核，
  以 COQ8A UMI / total RNA UMI 定義上下 25% 或上下 1/3；保留 cutoff ties。
  cutoff 相同時不強行劃組。
- Detected vs zero: 作 dropout-sensitive 對照，不等同真實生物學的完全無表現。
- Matching: none、0.10、0.20、0.30；在 sample×celltype 內，匹配
  log1p(total RNA counts) 與 log1p(detected ATAC peaks)，兩維差值均不得超過
  caliper。貪婪 nearest-neighbor matching，不放回，最多找 100 個候選鄰居。
  配對改變實際納入的核，因此它也改變估計所代表的細胞群。

舊流程的 anchor 為 author QC、≥3 vs 1、caliper 0.10、Pre，按細胞類型分開。
MuSC 只剩 3 對／1 donor，無 donor-level contrast；Fast/Slow/Intermediate
的 CAV3 FC 分別 0.978、0.934、0.822，p 分別 0.816、0.487、0.575。

### Count model

先在 donor×time×celltype×COQ8A group 內加總原始 peak counts。
**用全部 144,663 作者 peaks 計算 TMM library factors**，不是用選出的
候選 peaks 當作全 ATAC library size。候選過濾要求 total counts≥10、
至少兩位 donor 有計數（合併 high/low 後判定）。這是 outcome-blind
count-support filter；每組實際可檢驗的 peak 數另存表格。

edgeR 4.2.2 robust dispersion / quasi-likelihood model：

- matched: `~ donor + COQ8A_group`
- unmatched: `~ donor + within_donor_centered_mean_log_RNA_depth + COQ8A_group`

需要至少兩位有 high/low 的 donor、full-rank design 與 residual degrees
of freedom。無法估計的設定單獨記錄。既有 count-supported candidate peaks
用於 dispersion borrowing；TMM normalization 來源是全 ATAC。
RNA 平行使用 raw RNA library-size offset 的 donor QL model；RNA 結果
不是靠 ATAC 的 p 值判定。COQ8A/COQ9 不納入 myogenesis RNA FDR family。

ATAC p 在每個設定內分別按上述四個固定外部 scopes 做 BH；不是拿單一
peak 的 p 直接叫 q。**沒有再對 448 個 sensitivity settings 做一次 BH**；
最佳設定的 p/q 因而是探索中的 conditional statistics，不是搜尋程序的
整體錯誤率。每個 scope 的實際可檢驗數與全部結果均保留。

### Additional checks

開放比例為 peak raw count>0 的核比例。分別報 pooled detection FC、
equal-donor detection FC、donor paired t、精確 donor sign-flip；matched-nucleus
binomial p 另存，沒有代替 donor 主模型。

連續模型在 donor 內估 COQ8A normalized RNA 對 binary accessibility 的
partial correlation，控制 RNA depth、detected ATAC、mitochondrial fraction；
另加固定 myogenesis expression score。peak–RNA links 再控制 COQ8A，
並測試 target-excluded programme score。donor Fisher-z 等權合併、leave-one-out，
與最低偵測事件 3 改 1 的敏感度均保存。
補充 precision-aware random-effects Fisher-z meta-analysis 採 DL tau² 與
modified Hartung–Knapp variance floor；小樣本的 variance approximation
只作穩定性診斷，不用它取代另一個結果來挑顯著。

選定 CAV3 與另一 locus 另測試 donor-only / donor+RNA、exercise/rest、
leave-one-donor-out。這些都屬本次探索後的壓力測試。

## 3. CAV3：最有利條件與穩定程度

| Config | Population and definition | N donors; high/low nuclei | Count FC | p | Fusion-promoter q |
|---|---|---|---:|---:|---:|
| C296 | Slow Post, author QC, positive Q25, unmatched + RNA covariate | 6; 454/454 | 3.194 | 0.000860 | 0.01376 / 16 tested peaks |
| C324 | Same, depth cap95 | 6; 401/401 | 3.628 | 0.001471 | 0.02206 / 15 |
| C297 | C296 definition, depth match0.10 | 6; 176/176 | 1.329 | 0.08487 | 0.36776 / 13 |
| C298 | Same, match0.20 | 6; 257/257 | 1.051 | 0.72987 | 0.80354 / 14 |
| C299 | Same, match0.30 | 6; 318/318 | 1.058 | 0.66506 | 0.76737 / 15 |

C296 不加 RNA covariate 時 FC=1.173、p=0.207；C324 不加該 covariate 時
FC=0.956、p=0.747。這揭露的是**模型依賴**，不是擅自判定任一模型計算錯誤。
原始 Q25 組別 high 的平均 RNA depth 在各 donor 普遍較低；條件化於
RNA depth 後係數變大。C296 leave-one-donor-out 的 covariate model FC
仍為 2.83–3.50，但換深度處理方式，倍率大幅改變。

C297 的 detection 結果最容易直接閱讀：

| Donor | Paired nuclei per group | High open | Low open |
|---|---:|---:|---:|
| E | 50 | 20 | 18 |
| G | 33 | 12 | 7 |
| I | 42 | 7 | 5 |
| J | 37 | 14 | 13 |
| L | 5 | 3 | 1 |
| N | 9 | 6 | 3 |

合計 high 62/176=35.23%，low 47/176=26.70%，**pooled detection FC=1.319**。
六 donor 全同向，精確 sign-flip p=0.03125；23 mapped fusion promoters
的 detection-family q=0.71875。equal-donor FC=1.634，受到核數少的
L/N 比例影響較大，不能將它與 pooled FC 混寫。

相同 C297 RNA count model：CAV3 FC=1.061、p=0.790、q=0.972。
Slow Post 的 adjusted peak–CAV3 RNA link r≈0.044、p≈0.194。
所以這裡的支持是「某些條件下 CAV3 鄰近 accessibility 同向」，
尚未得到「同一條件 CAV3 RNA 也穩定提高」的完整鏈。

另一個已分析的 C2C12 分化參考 GSE224489，Cav3 在48h對growth medium
RNA FC=16.833、q=1.17×10⁻¹⁴，支持它與分化狀態相伴增加。
這是保留 CAV3 生物候選的另一理由；該對比本身並未操作 COQ8A。
細節與程式見 `CANDIDATE_EXPANSION_AND_EXTERNAL_DATA.md`。

整個 448-grid 中 CAV3 有 9 組 count FC>1 且 raw p<0.05，全部未配對；
5 組在 fusion-promoter family q<0.1。**所有可估的 matched count models
均未同時滿足 CAV3 FC>1、p<0.05。**這不否定前述 sign-flip，因為二者
檢驗的是不同的資料摘要與統計模型。

### MuSC continuous result

最低 3 個 COQ8A 偵測事件的分析保留 5 位 baseline donors，
5/5 partial-r>0，equal-donor r=0.0545、p=0.0196，範圍 q=0.938。
放寬到 1 event 納入第六位 L（12 顆核、2 顆 COQ8A+），其 r=-0.203，
合併 r=0.0112、p=0.813。依精確度加權的 6-donor 結果 r=0.0442、p=0.232。
因此不能把只保留 5 位的顯著 p 當成穩健的 MuSC 重現。

## 4. 全範圍與另一條候選

未匹配 Fast Post ≥3 vs1 的 HSMM module detection FC=1.140、p=0.00309，
但 high 組也較深。Fast Post Q33、caliper0.10 的全候選 module FC=1.0147、
p=0.00125，深度已相近。這說明平均開放幅度小不妨礙少数個別 peak
較大，亦說明未處理深度的整體倍率不能直接當作生物效應。

完整 candidate count screen 找到 **chr16:1311478–1312392**：

- Fast Post、COQ8A≥2 vs1、depth-cap95、caliper0.20：6 donors、490 對，
  FC=3.007、p≈8.19×10⁻⁶、q=0.0224（2,737 個可檢驗候選）。
- 同 QC、caliper0.10：FC=3.010、q=0.0271；0.30：FC=2.797、q=0.0916。
- 不加 cap95，配對的 FC 仍約 2.09–2.29，但全範圍 q=0.657–0.904。
- leave-one-donor-out 保持上升，FC=2.68–3.71，raw p 全部<0.003；
  每次重新檢驗全範圍後 q 並非均<0.1。
- 在已配對模型再加殘餘 RNA-depth 差異共變項後 p=0.254；
  兩種結果均保留，不依顯著與否挑選共變項。另存的 design condition
  number 受共變項尺度影響，不能單憑數值較高判定模型不可用。

原始候選註記將該位點分配為 CACNA1H 附近，到 CACNA1H 註記 TSS 距離約99 kb。
**它不是先前下降的 CACNA1H promoter peak。**此資料的 CACNA1H RNA
沒有通過相應計數/偵測支持條件，沒有可用的 cis RNA-link 結果。
補查不限制221基因的 GENCODE v48 後，此 peak 位於 UBE2I gene body，
並覆蓋 UBE2I-207 / ENST00000406620.6 的 TSS（0-based 1311870；
距 peak midpoint 65 bp）。這是更直接的 positional annotation，但仍非
已驗證 peak–RNA link。原來的99-kb是「到候選集合內CACNA1H的距離」，
**不是到全基因組最近gene的距離**。ATAC FC/p/q與檢驗範圍不變；
基因歸屬解釋以此修正為準。輸出在 `results/candidate_identity/`，
可用 `python scripts/72_candidate_identity_audit.py --gtf GENCODE_GTF_GZ` 重現。

### Fusion scope does not replace the myogenesis universe

12 genes 為 MSigDB 的 `GOBP_MYOBLAST_FUSION`（GO:0007520）與221全集交集，
不是完整 GO 集合只有12，也不是依ATAC顯著性選出的12。MYOD1/MYOG在此交集；
MYF5在 differentiation/positive-regulation scopes。MYOD1的一根500-bp
promoter未通過外部muscle-selectivity；MYOG共同候選peaks最近midpoint距TSS
3,758 bp，故沒有進500-bp promoter子集。三者仍在完整候選分析中。
其GSE208248共同peaks分別19/25/17；映射新GSE240061後分別11/12/11。
在新chr16 peak的C198設定，這三者鄰近可檢驗peaks沒有q<0.1。

## 5. 第二階段評級

| Route/claim | Grade | Decision |
|---|---|---|
| 只選最後剩 CAV3 一個 peak 的負調控集合，把 q=p 當強證據 | 1：不接受 | 功能方向不單一，且範圍過度依賴欲保留的候選 |
| 舊 fusion-promoter+external selectivity 路線，完整呈現敏感度 | 2：可作探索 | 有外部依據、有同規則其他候選；但門檻與窄範圍敏感、N=2來源 |
| 新 Slow Post 3.19–3.63 倍當作穩健獨立重現 | 1：不接受這項主張 | 效果依赖 RNA-depth 模型，匹配與 RNA 串接未穩定支持 |
| 新 Slow Post 配對六 donor 同方向、約1.32倍 | 2：可作探索 | 有同方向訊號，尚不具跨設定 FDR/轉錄連結穩定性；且是成熟肌核 |
| 固定外部範圍、完整 grid、donor模型、深度與來源敏感度的報告方式 | 3：方法可辯護 | 結論必須反映實際穩定程度，不保證 reviewer 接受有利的生物主張 |
| chr16 distal interval 的新發現 | 2：可作候選 | 有較強配對 ATAC 效果與範圍內 q，但 QC與gene linkage待解 |

本輪暫定：**CAV3 保留；不把新資料包裝成強力重現。**若只依這批公開資料
安排下一步，CAV3 適合作低成本 expression/targeted accessibility 候選，
不建議單靠 3 倍模型結果直接擴成大型 ChIP 機制計畫。
高低組條件一旦選定，global/module、single peak、RNA 與圖表都使用同一設定；
不同條件明列為 sensitivity，不再在下游偷偷換門檻。

## Figures and reading guide

**Figure 1 — External functional scope and GSE208248 CAV3 prioritization.**
Panel A includes every one of the 21 distinct fusion-associated promoter peaks
within the 500-bp midpoint-to-TSS window. Blue indicates external annotation
presence, not accessibility intensity. MYOD1 binding is annotation only, not an
additional selection gate. Panel B gives exact set sizes. Panel C shows pooled
detection fractions in 212 matched pairs from four libraries/two cell sources,
TSS≥2, COQ8A≥3 versus 1 UMI. Two-sided paired-nucleus binomial p and BH q across
the two retained peaks are shown. This is a retrospective prioritization scope.

**Figure 2 — CAV3 sensitivity across definitions and cell states.** Author-QC
settings only; each cell is normalized ATAC-count FC, with color showing log2FC.
The color scale saturates beyond ±2. Asterisks mark nominal donor QL p<0.05,
not FDR significance. Grey means no estimable CAV3 test. Left: all eligible
nuclei and RNA-depth covariate; right: caliper0.10 matched nuclei. Q25/Q33 refer
to tails among COQ8A-positive nuclei. Per-cell donor/sample sizes and q values
are provided in the complete source table. Panels target somewhat different
nucleus populations and are not interchangeable estimands.

**Figure 3 — Slow Post CAV3 follow-up.** Panel A connects each donor's low/high
peak-detection fraction under C297; labels give matched pairs. Donor E/G/I/J
are exercise, L/N rest. The displayed FC is pooled detection 62/47 and the p
is the two-sided exact donor sign-flip test. Detection q across the 23 mapped
fusion-promoter peaks is 0.71875. Panel B shows donor QL count-model FC and p
for C296–C299. Count-model q is calculated over 16/13/14/15 count-supported
fusion-promoter peaks, respectively; values are in the results table above.
Panel B point estimates summarize model sensitivity, not independent replicates.

**Supplement — Baseline MuSC donor sensitivity.** Technical-covariate partial
correlations for all six eligible donors under the minimum-one-event rule.
Each donor's nucleus and COQ8A-detected counts are shown; there is no exclusion
based on the sign of the correlation.

**Supplement — Distal interval sensitivity.** Panel A shows TMM-normalized
counts for chr16:1311478–1312392 under C198. The displayed p is donor QL and q
is BH over all 2,737 count-supported myogenesis-neighborhood peaks. Panel B
compares author-QC and depth-cap95 matched analyses, calipers0.10/0.20/0.30.
Each point's q uses its complete count-supported candidate family. This peak's
proximity to CACNA1H is an annotation, not a validated target assignment.

PNG and editable vector PDF: `figures/replication/`.
Figure data and numerical checks: `results/replication_summary/`.

## Reproduction

Obtain the author object and perform the inspection/exports described in
`GSE240061_FEASIBILITY_RESULTS.md` (scripts61–63). In the following, `WORK` is the
directory containing those exported matrices and metadata. `RDS` is the
decompressed author object. The large input and pseudobulk binaries remain
outside the repository. Python requirements are in `requirements.txt`;
R4.4.2, edgeR4.2.2, Matrix1.7-1, SeuratObject5.1.0 were used (full sessions saved).

```bash
python scripts/65_cav3_external_scope_grid.py
Rscript scripts/64_export_replication_background.R RDS WORK
python scripts/66_replication_sensitivity_inputs.py --work WORK --out results/replication_sensitivity
Rscript scripts/67_replication_donor_models.R . WORK results/replication_sensitivity
python scripts/68_replication_accessibility_effects.py --work WORK --out results/replication_sensitivity
python scripts/69_replication_continuous_links.py --work WORK --out results/replication_links --minimum-events 3
python scripts/69_replication_continuous_links.py --work WORK --out results/replication_links_relaxed --minimum-events 1
Rscript scripts/70_replication_model_stress_tests.R . WORK results/replication_stress
python scripts/71_summarize_replication.py
```

The last command alone regenerates the new figures and compact summaries from
committed result tables, with BH/FC consistency checks. It does not download or
rerun the large input pipeline. The complete 448-grid results are committed,
including unfavorable effects, failed models and low-information configurations.
Private C2C12 passage/time data remain outside this public repository.

## Method and biological references

- [GSE240061 integrated muscle multiome study](https://pmc.ncbi.nlm.nih.gov/articles/PMC12212352/).
- [Replicate-aware single-cell accessibility benchmark](https://www.nature.com/articles/s41467-024-53089-5).
- [edgeR user guide](https://www.bioconductor.org/packages/devel/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf).
- [CAV3/integrin involvement in myoblast fusion](https://pmc.ncbi.nlm.nih.gov/articles/PMC2710835/).

Functional association motivates testing CAV3; it does not mean every increase
in CAV3 dosage necessarily improves fusion. The public analysis is of COQ8A
expression association, not a COQ8A-overexpression intervention.

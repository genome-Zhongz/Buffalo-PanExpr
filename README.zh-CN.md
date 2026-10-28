# nonref-rna：非参考区域转录候选自动分析

输入双端 RNA-seq、样本表以及已有的水牛扩展参考文件包，一条命令运行比对、组装、定量、候选筛选和序列导出。可选启用同源注释、差异分析和透明背景图。

## 安装与运行

```bash
conda env create -f environment.yml
conda activate nonref-rna
python -m pip install .

nonref-rna init --reference-bundle /path/to/buffalo_annotation --out config.local.yaml
# 修改配置，参考 examples/samples.tsv 填写实际双端文件。
nonref-rna validate --config config.local.yaml --samples samples.tsv
nonref-rna run --config config.local.yaml --samples samples.tsv --out results/run1 --dry-run
nonref-rna run --config config.local.yaml --samples samples.tsv --out results/run1
```

完整运行要求 Linux。Python 主程序需要 PyYAML；其他生信工具由环境文件安装。环境文件不是已经求解的版本锁文件，实际版本写入运行日志。

## 输入文件

样本表为真实 TAB 分隔，必需列：`sample、r1、r2`。可加 `tissue、breed、group、individual、strandedness`。FASTQ 支持未压缩和 gzip；路径相对于样本表解析。样本和文件不能重复，不自动排除任何品种或组织。

`reads.cleaned: true` 复用质控数据；`false` 自动运行 fastp。链特异性支持 auto/rf/fr/unstranded。

当前版本复用预先建立的 `buffalo_annotation`：extended.fa、guide.gtf、nonref.guide.gtf 和 reference.normalized.gtf。具体布局和参数覆盖见英文 README 与 examples/config.yaml。

## 主要结果

- `candidates/nonreference_candidates.tsv`：非参考候选、表达和结构证据。
- `candidates/longest_transcript_per_gene.fa`：每个候选基因一条最长转录本。
- `candidates/all_transcripts.fa`：保留所有转录本。
- `candidates/gene_transcript_mapping.tsv`：基因、转录本、坐标与长度对应。
- `matrices/gene_TPM.tsv`、`gene_counts.raw.tsv`：表达矩阵和整数计数。
- `homology/gene_homology.tsv`：可选的本地 DIAMOND BLASTX 注释。
- `report/report.html`、`summary.json`：可追溯汇总。
- `figures/`、`differential_expression/`：启用后生成的 R 图和差异结果。

默认候选：RNA 实际组装位点位于非参考 contig，且至少一个样本 TPM ≥1。可提高 `expression.min_samples`。组织表达重叠另外按“组织平均 TPM >1”统计，跨组织去重；各组织数量直接相加不是唯一基因总数。

RNA 组装、表达、AUGUSTUS 结构匹配与蛋白同源证据分列呈现。不能把非参考位置、gffcompare 的 u 类或一次蛋白匹配直接解释成新功能基因。FASTA 是剪接后的核酸转录本，不是蛋白质，也不保证有完整 CDS。

## 可选功能

同源注释：设置 `homology.enabled: true` 和本地 `database`。保留蛋白描述、登录号、一致性、覆盖度和 E-value。“高一致性且覆盖较广”要求一致性 ≥80%、目标蛋白覆盖度 ≥80%、比对长度 ≥50 aa、E-value ≤1e-10；这不是转录本完整度指标。只有数据库标题带 GN 字段时提取符号，否则报告蛋白描述和登录号。

`homology/homology_summary.md` 自动统计当前表达候选的蛋白匹配数量；`primary_Pvalue/homology_summary.md` 单独统计原始 P 值差异候选。两者均提供中英文描述和 JSON，HTML 报告也显示汇总。例如，实际结果为 48/42/25 时会输出：“48 个候选中，42 个检出蛋白同源匹配，其中 25 个达到较高一致性和较广蛋白覆盖标准。”数量按基因统计，不按匹配蛋白或转录本条数累计；未注释和搜索后未检出分开统计。

差异分析：设置 `analysis.differential: true`、control/treatment 和样本 group。使用原始 counts、CPM 过滤、TMM 和 edgeR robust QL 检验。当前只支持非配对两组比较，每组至少两个独立生物学样本；重复测量、配对设计需另行建模。原始 P 和 FDR 结果分别输出。

**主结果采用原始 P<0.05 且 |log2FC|>1。优先查看 `primary_Pvalue/`：**

- `nonreference_candidates_P05.csv`：按 P 值排序的非参考候选，含方向、可用的同源蛋白描述与转录本 ID。
- `nonreference_candidates_P05.longest.fa`：这些候选每基因一条最长转录本；另有全部转录本与 GTF。
- `Volcano_Pvalue.png`、`MDS_TPM.png`：启用绘图后生成的透明主图；候选不足时不虚构 MDS。
- `all_genes_P05.csv`：含参考和非参考的全部原始 P 候选。

FDR 结果只在完整 `differential_expression/` 目录中保留供参考。差异候选不再用 TPM 阈值二次过滤，以免丢失通过原始 counts 检验的基因。

绘图：设置 `analysis.plots: true`，R 输出透明 PNG/PDF，包括表达数量、组织交集、TPM PCA；有差异分析时输出火山图和候选 TPM MDS。筛选后 MDS 仅作描述性展示。

## 断点续跑与已有结果

重复同一运行命令即可续跑。改变线程等资源设置不改变生物学输入；改变样本、参考、工具或核心参数须使用新目录，防止混用。文件默认按路径、大小、修改时间核验，代码用 SHA256；不对数百 GB 的 FASTQ 默认全量哈希。

```bash
nonref-rna summarize --work /path/to/completed/bam-new --out results/imported
nonref-rna status --out results/run1
```

历史目录必须保留 samples.tsv、表达与元数据矩阵、最终计数 GTF、非参考合并 GTF 和其 FASTA、各样本 quant/gene_abundance.tsv，以及源运行配置、输入清单和原始扩展参考。导入核对元数据坐标、逐样本 TPM 来源，并从参考外显子重建序列逐条比对；不能把不同组装的 BFRNA 编号直接混合。启用差异分析还需要原始 counts。

FASTQ 预检检查首对记录，自动链特异性推断额外检查 pilot 部分配对；不会另做全文件 read-ID 配对扫描。请保证输入文件配对正确，并检查质控与比对日志。

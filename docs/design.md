# nonref-rna：Linux 命令行软件设计

状态：用户已于 2026-10-28 确认设计；实现与验证记录见 VALIDATION.md。

## 目标

将已用于水牛转录组分析的流程整理成可上传 GitHub 的 Python 命令行项目。输入 RNA-seq 双端 FASTQ、样本表和已构建的参考文件包，自动输出非参考区域的 RNA 组装模型、表达证据、转录本 FASTA 和可选蛋白同源注释。输出称为“非参考候选基因位点”，不自动宣称发现全新功能基因。

采用 Linux 命令行、配置文件和样本表的形式，已由用户选定。Python 负责任务调度、数据校验和结果汇总；R 负责可选差异分析与静态科研绘图。外部程序保持 fastp、HISAT2、SAMtools、StringTie2、featureCounts、gffcompare、gffread、DIAMOND。

## 现有代码中可复用的部分

- HISAT2 --dta 比对、SAMtools 排序、CRAM 保存。
- 参考外显子重叠法推断链特异性，并同步设置 HISAT2、StringTie 和 featureCounts 参数。
- 每样本 StringTie 组装，合并实际观察到的模型，提取非参考 contig 上的 RNA 模型。
- 原始参考注释与 RNA 非参考模型组成定量 GTF；StringTie -e 定量，featureCounts 输出整数 counts。
- TPM/counts 矩阵、样本质量报告、gffcompare 对 AUGUSTUS 的结构支持。
- 断点续跑、多样本并行、磁盘余量控制、文件与参数溯源。
- 候选基因全部转录本和每基因最长转录本导出；DIAMOND BLASTX 同源注释。

通用化时去掉固定的组织名称、样本目录发现规则和服务器程序路径。资源修改与参考修改的续跑逻辑分开；不同参考、样本或生物学参数不得沿用旧结果。

## 输入与运行方式

样本表必需列：sample、r1、r2。可选列：tissue、breed、group、individual、strandedness。样本 ID 必须唯一；路径相对于样本表解析。首版支持双端 RNA-seq，暂不支持单细胞、长读长或单端数据。

FASTQ 可为原始或已质控数据，由配置显式标记；原始数据自动运行 fastp，已质控数据复用。

参考输入为当前构建的文件包：

- rnaseq_reference/extended.fa：扩展参考；
- rnaseq_reference/guide.gtf：参考和预测模型的组合引导注释；
- reference/reference.normalized.gtf：原始参考注释；
- rnaseq_reference/nonref.guide.gtf：已有非参考预测模型；
- 非参考 contig 清单或前缀；
- 可选 combined_db.dmnd 与现有 AUGUSTUS/BLASTP 证据。

初始化命令根据文件包生成配置并解析 PATH 中的软件，不要求其他用户复制本服务器软件安装路径。公开示例使用占位路径；个人服务器配置留在 Git 忽略的本地目录。数据库、FASTQ、BAM、CRAM 和实际参考 FASTA 不加入 GitHub 仓库。

预期命令接口：

```bash
nonref-rna init --reference-bundle /path/to/buffalo_annotation --out config.yaml
nonref-rna validate --config config.yaml --samples samples.tsv
nonref-rna run --config config.yaml --samples samples.tsv --out results
```

同一命令重新运行可续跑。提供 --dry-run 检查参数、软件和输入并打印流程，不启动全量比对。另提供从本软件已完成运行目录重新汇总候选的命令；对历史项目仅在注释、参考和定量来源验证通过后导入，不能只按 BFRNA 名称盲目拼接。

## 自动流程

1. 校验工具、参考/GTF contig 一致性、样本 ID、FASTQ 配对及资源配置。
2. 按配置执行或跳过 fastp；记录质控来源。
3. 构建或复用经过核验的扩展参考 HISAT2 索引；可选比对原始参考作为比较。
4. 判断或读取链特异性，逐样本比对及组装。
5. 合并观察到的 RNA 模型，明确提取非参考 contig 上的模型，生成统一定量注释。
6. 固定模型 TPM 定量及整数 read/fragment counts 汇总。
7. 按可配置表达规则提取非参考候选；输出按样本、组织/分组平均和跨样本去重的不同统计口径。
8. 与已有 AUGUSTUS 预测进行结构比较；可选 DIAMOND BLASTX 及蛋白来源注释。
9. 导出 GTF、FASTA、映射表、证据表、汇总报告；有分组时可选运行差异分析和图形输出。

## 候选标准与解释边界

默认表达标准沿用现有证据表：TPM >=1，至少一个样本满足；最少支持样本数和阈值可调。组织 UpSet 独立使用“组织平均 TPM >1”，报告中明确区分 > 与 >=，不混用组织和样本口径。

证据分列报告：RNA 组装、表达支持、AUGUSTUS gffcompare 类别、同源匹配一致性、查询覆盖、蛋白覆盖、比对长度及 E-value。可筛选组合证据，但不将任何蛋白命中或 gffcompare 的 u 类直接解释为新基因。

将与已知蛋白高度相似、多个候选对应同一蛋白、部分匹配、含 N 的序列和无法确定编码潜能等情况显式列出；不声称仅凭这些证据已经确认假基因或新功能。

同源库与 NCBI 在线 ID 查询均为可选项；默认不上传用户序列到第三方网站。未配置数据库时，蛋白注释状态为“未运行”，不写成“无同源”。

差异分析使用原始 counts 和 edgeR，表达过滤使用 CPM，TMM 标准化。TPM 用于表达描述与 log2(TPM+1) 的可选 PCA/MDS。FDR 与原始 P 筛选分开输出，不将原始 P 候选标成 FDR 显著。基于已筛选候选的 MDS 明确标注为描述性展示。

## 输出

```text
results/
  run_manifest.json
  status.json
  logs/
  annotation/
  matrices/gene_TPM.tsv
  matrices/gene_counts.raw.tsv
  candidates/nonreference_candidates.tsv
  candidates/all_transcripts.fa
  candidates/longest_transcript_per_gene.fa
  candidates/gene_transcript_mapping.tsv
  candidates/selected_models.gtf
  homology/                       # 配置数据库时生成
  differential_expression/       # 用户启用且分组有效时生成
  figures/                       # 可选 R 绘图
  report/summary.json
  report/report.html
```

图形可选包括表达数量、组织重叠和 PCA/MDS；有差异分析时增加火山图。所有数字可追溯到源表，默认支持透明导出。

## 资源、运行安全和复现

总体 CPU 预算可配置，并发样本数与每阶段线程数受总预算约束；记录磁盘需求和不足提示。程序仅清理自身工作目录的临时文件。失败记录具体阶段和命令，重跑先核验输入与产物，禁止悄悄跨参考复用。

首版交付源码目录与 GitHub-ready ZIP，包含 pyproject.toml、环境文件、英文/中文 README、示例配置和样本表、测试、GitHub Actions 和引用信息模板。软件许可证由作者发布时选择；不擅自替用户对研究数据或第三方数据库重新授权。

## 验收

- 单元测试覆盖非参考筛选、去重统计、空结果、基因/转录本映射、链特异性和断点续跑失效规则。
- 小型合成数据验证外部程序串联及候选输出；测试数据明确标记为合成，不能当作生物学发现。
- 利用现有水牛结果复核表达统计及48候选序列提取/同源注释的兼容性，不重新运行全部73个文库。
- 报告真实测试范围；未执行的全量测试不能宣称通过。
- 交付打包文件，不自动发布 GitHub 仓库或上传用户数据。

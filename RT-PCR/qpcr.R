#' 自动化处理qPCR数据（自定义对照组，多Test组，3个生物学重复，内参校正）
#'
#' @param qpcr_df data.frame，qPCR原始数据，列要求：
#'                - Gene：基因名（含内参基因，如GAPDH/ACTB）
#'                - Group：分组名（任意命名，如"Control"、"DrugA"、"DrugB"）
#'                - Replicate：生物学重复（1/2/3）
#'                - Ct：原始Ct值（数值型，不可含NA/异常值）
#' @param ref_gene 字符型，内参基因名（如"GAPDH"）
#' @param control_group 字符型，自定义对照组名称（如"Control"、"NC"、"Mock"）
#' @param sig_level 数值型，统计检验显著性水平（默认0.05）
#' @param plot_show 逻辑型，是否展示可视化图（默认TRUE）
#' @param plot_save 逻辑型，是否保存可视化图（默认FALSE）
#' @param save_path 字符型，图片保存路径（默认当前目录）
#'
#' @return 列表，包含：
#'         1. raw_data：输入的原始数据（含质控标记）
#'         2. norm_data：归一化后的数据（ΔCt）
#'         3. rel_expr：相对表达量（2^-ΔΔCt）
#'         4. stats_result：统计检验结果
#'         5. plot：可视化ggplot对象（若plot_show=TRUE）
#'
#' @examples
#' # 模拟数据运行示例（自定义对照组为"Control"）
#' set.seed(123)
#' test_df <- data.frame(
#'   Gene = rep(c("GeneA", "GeneB", "GAPDH"), each = 12),  # 2个目标基因+1个内参
#'   Group = rep(c("Control", "DrugA", "DrugB", "DrugC"), each = 3, times = 3),  # 对照组为Control
#'   Replicate = rep(1:3, times = 12),
#'   Ct = c(
#'     rnorm(3, 25, 0.5), rnorm(3, 27, 0.5), rnorm(3, 26, 0.5), rnorm(3, 28, 0.5),  # GeneA
#'     rnorm(3, 24, 0.5), rnorm(3, 26, 0.5), rnorm(3, 25, 0.5), rnorm(3, 27, 0.5),  # GeneB
#'     rnorm(3, 20, 0.2), rnorm(3, 20.1, 0.2), rnorm(3, 20.2, 0.2), rnorm(3, 20.1, 0.2)  # GAPDH
#'   )
#' )
#' # 运行函数（指定对照组为"Control"）
#' result <- qPCR_analysis(qpcr_df = test_df, ref_gene = "GAPDH", control_group = "Control")
#'
#' @import dplyr
#' @import ggplot2
#' @import rstatix
#' @import ggpubr
#' @import RColorBrewer
#' @export
qPCR_analysis <- function(qpcr_df,
                          ref_gene,
                          control_group,  # 新增：自定义对照组参数
                          sig_level = 0.05,
                          plot_show = TRUE,
                          plot_save = FALSE,
                          save_path = getwd()) {
  
  # 加载必需包（无则自动安装）
  required_pkgs <- c("dplyr", "ggplot2", "rstatix", "ggpubr", "RColorBrewer")
  for (pkg in required_pkgs) {
    if (!require(pkg, character.only = TRUE)) {
      install.packages(pkg, dependencies = TRUE)
      library(pkg, character.only = TRUE)
    }
  }
  
  # ===================== 步骤1：数据质控（新增对照组校验） =====================
  cat("开始数据质控...\n")
  # 检查列是否完整
  required_cols <- c("Gene", "Group", "Replicate", "Ct")
  if (!all(required_cols %in% colnames(qpcr_df))) {
    stop(paste("数据框缺少必需列！必需列：", paste(required_cols, collapse = ", ")))
  }
  
  # 检查数据类型
  qpcr_df <- qpcr_df %>%
    mutate(
      Gene = as.character(Gene),
      Group = as.character(Group),
      Replicate = as.factor(Replicate),
      Ct = as.numeric(Ct)
    )
  
  # 新增：检查自定义对照组是否存在
  if (!control_group %in% unique(qpcr_df$Group)) {
    stop(paste("指定的对照组'", control_group, "'未在数据中找到！\n已有的分组：", 
               paste(unique(qpcr_df$Group), collapse = ", ")))
  }
  
  # 检查内参基因是否存在
  if (!ref_gene %in% unique(qpcr_df$Gene)) {
    stop(paste("内参基因'", ref_gene, "'未在数据中找到！"))
  }
  
  # 检查每组每个基因的重复数是否为3
  rep_check <- qpcr_df %>%
    group_by(Gene, Group) %>%
    summarise(n = n(), .groups = "drop") %>%
    filter(n != 3)
  if (nrow(rep_check) > 0) {
    warning(paste("以下基因-分组的重复数不是3个：\n",
                  paste(rep_check$Gene, rep_check$Group, paste("n=", rep_check$n), sep = " | ", collapse = "\n")))
  }
  
  # 检查Ct值异常（Ct<10或Ct>40视为异常）
  qpcr_df <- qpcr_df %>%
    mutate(Ct_quality = ifelse(Ct < 10 | Ct > 40, "异常", "正常"))
  abnormal_ct <- qpcr_df %>% filter(Ct_quality == "异常")
  if (nrow(abnormal_ct) > 0) {
    warning(paste("检测到异常Ct值（<10或>40）：\n",
                  paste(abnormal_ct$Gene, abnormal_ct$Group, paste("Ct=", abnormal_ct$Ct), sep = " | ", collapse = "\n")))
  }
  
  # ===================== 步骤2：计算ΔCt（目标基因Ct - 内参基因Ct） =====================
  cat("开始计算ΔCt...\n")
  # 提取内参基因的Ct值（按分组+重复匹配）
  ref_ct <- qpcr_df %>%
    filter(Gene == ref_gene) %>%
    select(Group, Replicate, ref_Ct = Ct)
  
  # 合并内参Ct值，计算ΔCt
  norm_data <- qpcr_df %>%
    filter(Gene != ref_gene) %>%  # 排除内参基因，仅处理目标基因
    left_join(ref_ct, by = c("Group", "Replicate")) %>%
    mutate(DeltaCt = Ct - ref_Ct) %>%  # ΔCt = 目标基因Ct - 内参Ct
    select(Gene, Group, Replicate, Ct, ref_Ct, DeltaCt, Ct_quality)
  
  # ===================== 步骤3：计算ΔΔCt和相对表达量（2^-ΔΔCt）（适配自定义对照组） =====================
  cat("开始计算相对表达量（2^-ΔΔCt）...\n")
  # 计算自定义对照组每个基因的平均ΔCt
  ctrl_mean_dct <- norm_data %>%
    filter(Group == control_group) %>%
    group_by(Gene) %>%
    summarise(Ctrl_mean_DeltaCt = mean(DeltaCt, na.rm = TRUE), .groups = "drop")
  
  # 计算ΔΔCt和相对表达量（实验组ΔCt - 对照组平均ΔCt）
  rel_expr <- norm_data %>%
    left_join(ctrl_mean_dct, by = "Gene") %>%
    mutate(
      DeltaDeltaCt = DeltaCt - Ctrl_mean_DeltaCt,  # 适配自定义对照组
      RelExpr = 2^(-DeltaDeltaCt)  # 2^-ΔΔCt 相对表达量
    ) %>%
    select(Gene, Group, Replicate, Ct, ref_Ct, DeltaCt, DeltaDeltaCt, RelExpr, Ct_quality)
  
  # ===================== 步骤4：统计检验（实验组vs 自定义对照组） =====================
  cat("开始统计检验...\n")
  # 定义实验组（排除对照组）
  test_groups <- setdiff(unique(rel_expr$Group), control_group)
  stats_result <- data.frame()
  
  if (length(test_groups) == 1) {
    # 单实验组：与自定义对照组做t检验
    stats_result <- rel_expr %>%
      group_by(Gene) %>%
      t_test(RelExpr ~ Group, ref.group = control_group, p.adjust.method = "BH") %>%
      add_significance("p.adj")
  } else {
    # 多实验组：ANOVA + Tukey HSD事后检验（仅对比自定义对照组）
    anova_res <- rel_expr %>%
      group_by(Gene) %>%
      anova_test(RelExpr ~ Group) %>%
      add_significance("p.adj")
    
    tukey_res <- rel_expr %>%
      group_by(Gene) %>%
      tukey_hsd(RelExpr ~ Group) %>%
      filter(str_detect(contrast, control_group))  # 仅保留与自定义对照组的对比
    
    stats_result <- left_join(anova_res, tukey_res, by = "Gene")
  }
  
  # ===================== 步骤5：可视化（适配自定义对照组） =====================
  cat("开始绘制可视化图...\n")
  # 整理绘图数据
  plot_data <- rel_expr %>%
    filter(Ct_quality == "正常") %>%  # 排除异常Ct值
    mutate(Group = factor(Group, levels = c(control_group, test_groups)))  # 固定分组顺序（对照组在前）
  
  # 绘制箱线图（每个基因单独面板）
  p <- ggplot(plot_data, aes(x = Group, y = RelExpr, fill = Group)) +
    geom_boxplot(width = 0.6, outlier.shape = 21, outlier.size = 2) +  # 箱线图
    geom_jitter(width = 0.1, size = 1.5, alpha = 0.6) +  # 散点展示重复
    facet_wrap(~Gene, scales = "free_y") +  # 按基因分面
    # 配色：对照组为灰色，实验组为彩色
    scale_fill_manual(values = c(
      control_group = "gray80", 
      setNames(brewer.pal(min(length(test_groups), 8), "Set1"), test_groups)  # 适配多实验组配色
    )) +
    # 显著性标记（对比自定义对照组）
    stat_compare_means(
      ref.group = control_group,
      method = ifelse(length(test_groups) == 1, "t.test", "anova"),
      label = "p.adj.signif",
      hide.ns = TRUE,
      size = 3.5
    ) +
    labs(
      x = "分组",
      y = "相对表达量 (2^-ΔΔCt)",
      title = "qPCR相对表达量分析",
      subtitle = paste("内参基因：", ref_gene, " | 对照组：", control_group)
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 12),
      axis.text = element_text(size = 10),
      axis.title = element_text(size = 12),
      legend.position = "none",
      strip.text = element_text(size = 11, face = "bold")
    )
  
  # 展示图片
  if (plot_show) {
    print(p)
  }
  
  # 保存图片
  if (plot_save) {
    ggsave(
      filename = paste0(save_path, "/qPCR_analysis_", Sys.Date(), ".pdf"),
      plot = p,
      width = 8 + length(test_groups),
      height = 6,
      dpi = 300
    )
    cat(paste("图片已保存至：", save_path, "/qPCR_analysis_", Sys.Date(), ".pdf\n"))
  }
  
  # ===================== 返回结果 =====================
  result <- list(
    raw_data = qpcr_df,          # 原始数据（含质控）
    norm_data = norm_data,       # ΔCt归一化数据
    rel_expr = rel_expr,         # 相对表达量（2^-ΔΔCt）
    stats_result = stats_result, # 统计检验结果
    plot = p                     # 可视化图对象
  )
  
  cat(paste("qPCR数据处理完成！对照组为：", control_group, "\n"))
  return(result)
}

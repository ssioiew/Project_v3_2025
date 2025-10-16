# ------------------------------------------------------
# File: Project v.3.3 2025
# Purpose: Brief description of the script
# Author: Pirada Naewkam
# Date: 2025-10-09
#
# Description:
#   Detailed description or steps
#
# Input: bracken file
# Output: visualizations and statistical test
# ------------------------------------------------------


# Step 0: Setup Environment ----------------------------------------------

library(tidyverse)
library(purrr)
library(readxl)
library(forcats)
library(RColorBrewer)
library(ggplot2)
library(scales)
library(vegan)
library(ggsignif)
library(ggpubr)
library(ggrepel)
library(MicrobiotaProcess)
library(cowplot)
library(patchwork)
library(stringr)
library(rstatix)
library(multcompView)
library(showtext)
font_add_google("Lora", "lora_font")  # "lora_font" คือชื่อเล่นที่เราจะใช้เรียกใน R
font_add(family = "times_new_roman", # ตั้งชื่อเล่นที่เราจะใช้ใน R
         regular = "times.ttf",     # ชื่อไฟล์ฟอนต์ปกติ
         bold = "timesbd.ttf",      # ชื่อไฟล์ฟอนต์ตัวหนา
         italic = "timesi.ttf",     # ชื่อไฟล์ฟอนต์ตัวเอียง
         bolditalic = "timesbi.ttf") # ชื่อไฟล์ฟอนต์หนาและเอียง
showtext_auto()

# File Paths
DATA_FOLDER   <- "C:/Old Volume/Work/Project Work/Project_v3_2025/Kraken Data/trimmomatic2/Map/"
METADATA_FILE <- "C:/Old Volume/Work/Project Work/Project_v3_2025/sample_metadata.csv"

# Plotting Dimensions & Settings
PUB_WIDTH <- 7
PUB_HEIGHT <- 4.5
PUB_DPI <- 300

# สร้าง function palette สีตาม Location
create_color_palette <- function(group_vector, light_palette = "Set2", dark_palette = "Dark2") {
  
  # สร้างลำดับมาตรฐาน
  group_order <- sort(unique(group_vector))
  num_groups <- length(group_order)
  
  # สร้าง palette สี
  light_colors <- brewer.pal(n = num_groups, name = light_palette)
  dark_colors <- brewer.pal(n = num_groups, name = dark_palette)
  
  # สร้าง name vector
  names(light_colors) <- group_order
  names(dark_colors) <- group_order
  
  # ส่งค่ากลับไปให้ object
  return(list(light = light_colors, dark = dark_colors))
}


# Step 1: Load and Process Raw Data -------------------------------------------------------

# โหลดข้อมูล bracken ทั้งหมด
retrieve_data <- list.files(paste0(DATA_FOLDER), pattern = "*.bracken.genus.txt", full.names = TRUE)

# ดึง sample id ออกมาจากชื่อไฟล์โดยตรง
sample_ids_from_files <- sub(".bracken.genus.txt$", "", basename(retrieve_data))

sample_ids_from_files <- sub("\\..*", "", basename(retrieve_data))

# โหลด data และผูก sample id เข้าไป
bracken_data <- lapply(seq_along(retrieve_data), function(i) {
  df <- read.delim(retrieve_data[i])
  df$sample_id <- sample_ids_from_files[i]  # ใช้ id จากไฟล์
  return(df)
}) %>% bind_rows()

# สร้าง metadata frame ที่สมบูรณ์
sample_info <- read_csv(METADATA_FILE) %>% 
  mutate(Location = case_when(
    
    # ใช้ str_detect() เพื่อมองหา "คำสำคัญ" ในชื่อแหล่งโบราณคดี
    str_detect(`Archaeological site`, "NRW") ~ "NRW",
    str_detect(`Archaeological site`, "TamOngba") ~ "TamOngba",
    str_detect(`Archaeological site`, "KK") ~ "KK",
    str_detect(`Archaeological site`, "PMN")  ~ "PMN",
    str_detect(`Archaeological site`, "BanKao") ~ "BanKao",
    str_detect(`Archaeological site`, "BanNaDi")  ~ "BanNaDi",
    TRUE ~ "Unknown")
  )

metadata_df <- data.frame(sample_id = sample_ids_from_files) %>% 
  left_join(sample_info, by = c("sample_id" = "Sample ID"))


# Step 2: Centralized Data Wrangling & Pre-computation --------------------

# คำนวณ Total Reads ของแต่ละ sample
sample_totals <- bracken_data %>% 
  group_by(sample_id) %>% 
  summarise(total_reads = sum(new_est_reads))

# Join total_reads กลับเข้าไป และคำนวณ % จริงของแต่ละ sample
bracken_data <- bracken_data %>% 
  left_join(sample_totals, by = "sample_id") %>% 
  mutate(percent_reads = (new_est_reads / total_reads) * 100) %>% 
  
  # จัดลำดับ taxa ตาม abundance รวมทั้งหมด
  mutate(name = fct_reorder(name, new_est_reads, .fun = sum, .desc = TRUE))

# สร้าง OTU Table (Wide Format) สำหรับ Diversity Analysis
bracken_data_wide <- bracken_data %>% 
  select(sample_id, name, new_est_reads) %>% 
  pivot_wider(names_from = name, values_from = new_est_reads, values_fill = list(new_est_reads = 0))

# แปลงเป็น matrix ให้ vegan
bracken_matrix <- bracken_data_wide %>% 
  select(-sample_id) %>% 
  as.matrix()
rownames(bracken_matrix) <- bracken_data_wide$sample_id

# คำนวณ Alpha Diversity และสร้าง dataframe
shannon_diversity <- diversity(bracken_matrix, index = "shannon")

# สร้าง dataframe สำหรับพล็อต
alpha_diversity_df <- data.frame(
  sample_id = names(shannon_diversity),
  shannon_index = shannon_diversity
) %>%
  left_join(metadata_df, by = "sample_id")  # join กับข้อมูล metadata (Location)

# เรียกใช้งาน function สี ตาม Location
location_palettes <- create_color_palette(alpha_diversity_df$Location)


# Step 3: Stacked barplot ---------------------------------------------------------

# Exploratory Step - Find the best "n"
# taxa_summary <- bracken_data %>%
#   group_by(name) %>%
#   summarise(total_abundance = sum(percent_reads)) %>%
#   arrange(desc(total_abundance)) %>%
#   mutate(
#     percent_of_total = total_abundance / sum(total_abundance) * 100,
#     cumulative_percent = cumsum(percent_of_total)
#   )
# print(taxa_summary, n = 25) # sweet spot is n <= 17
# 
# # Set choosen "n"
# CHOSEN_N <- 14
# 
# # สร้างข้อมูลสำหรับการพล็อต โดยจัดเรียงข้อมูลของแกน x ตาม Location
# stacked_data <- bracken_data %>% 
#   left_join(select(metadata_df, sample_id, Location), by = "sample_id") %>% 
#   
#   # เรียงลำดับ sample_id ใหม่ ตาม Location
#   mutate(sample_id = fct_reorder(sample_id, as.numeric(factor(Location)))) %>% 
#   
#   # สร้างข้อมูลสำหรับพล็อต: รวม Taxa ที่ไม่ใช่ Top 14 เป็น "Other"
#   mutate(name_lumped = fct_lump_n(name, n = CHOSEN_N, w = percent_reads, other_level = "Other")) %>% 
#   
#   # Aggregate data เพื่อให้ Other เป็นแถวเดียวต่อ sample
#   group_by(sample_id, Location, name_lumped) %>% 
#   summarise(total_percent = sum(percent_reads), .groups = "drop") %>% 
#   
#   # จัดเรียง legend ให้ถูกต้อง
#   mutate(
#     name_lumped = fct_reorder(name_lumped, total_percent, .fun = sum, .desc = TRUE),   # เรียงตาม abundance รวม
#     name_lumped = fct_relevel(name_lumped, "Other", after = Inf)   # ย้าย Other ไปไว้ท้ายสุด
#   )
# 
# taxa_levels <- levels(stacked_data$name_lumped)
# 
# # กำหนดสี special case
# manual_colors <- c(
#   "Homo"  = "#D8C3A5",
#   "Other" = "gray90"   
# )
# 
# # หาชื่อ taxa ที่เหลือ ที่จะใส่ palette สี
# remaining_taxa <- setdiff(taxa_levels, names(manual_colors))
# 
# # สร้าง Palette สีสำหรับ taxa ที่เหลือ
# color_generator <- c(RColorBrewer::brewer.pal(12, "Paired"), RColorBrewer::brewer.pal(12, "Set3"))
# palette_for_remaining <- color_generator[1:length(remaining_taxa)]
# 
# # ติดป้าย names ให้กับ palette สี taxa ที่ต้องการ
# names(palette_for_remaining) <- remaining_taxa
# 
# # รวมทุกสีเข้าด้วยกัน
# final_color_palette <- c(palette_for_remaining, manual_colors)
# 
# # พล็อตกราฟ
# stacked_barplot_v2 <-
#   ggplot(stacked_data, aes(x = sample_id, y = total_percent, fill = name_lumped)) +
#   geom_col(position = "stack", width = 0.75, color = "black", linewidth = 0.1) +
#   
#   facet_grid(~ Location, scales = "free_x", space = "free_x") +
#   
#   scale_y_continuous(expand = c(0, 0), limits = c(0, 101)) +  # ทำให้แท่งกราฟเริ่มที่ 0 พอดี และเผื่อที่ด้านบนเล็กน้อย
#   scale_fill_manual(
#     values = final_color_palette,
#     name = "Taxa",
#     breaks = levels(stacked_data$name_lumped)   # ทำให้ลำดับใน Legend ตรงกับลำดับใน Factor (Abundance สูงสุดอยู่บนสุด)
#   ) +
#   
#   labs(
#     title = "Relative Abundance of Dominant Taxa",
#     subtitle = "Top 14 most abundant genera are shown; others are grouped.",
#     x = "Sample ID",
#     y = "Relative Abundance"
#   ) +
#   
#   theme_bw(base_size = 22) + # เพิ่ม base_size เพื่อให้ตัวอักษรโดยรวมใหญ่ขึ้น
#   theme(
#     plot.title = element_text(face = "bold", hjust = 0.5),
#     plot.subtitle = element_text(hjust = 0.5, face = "italic"),
#     
#     axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), # ปรับมุม
#     axis.title = element_text(face = "bold"),
#     axis.line = element_line(linewidth = 0.4),
#     axis.ticks = element_line(linewidth = 0.4),
#     
#     legend.title = element_text(face = "bold"),
#     legend.background = element_blank(), # เอาพื้นหลัง legend ออก
#     legend.key.size = unit(0.6, "cm"),
#     
#     panel.grid.major = element_blank(), # เอาเส้นกริดแนวตั้งออกเพื่อให้แท่งบาร์เด่นขึ้น
#     panel.grid.minor = element_blank(),
#     panel.border = element_rect(color = "black", linewidth = 1),
#     panel.spacing.x = unit(0.1, "lines"),   # ลดช่องว่างของ Facet
#     
#     # ปรับแต่งแถบหัวเรื่อง Facet
#     strip.background = element_rect(fill = "gray85", color = "black"),
#     strip.text.x = element_text(face = "bold", size = 14),
#     
#     text = element_text(family = "times_new_roman")
#   )
# 
# print(stacked_barplot_v2)

# Export for publication
# output_folder <- "C:/Users/patch/Desktop/"
# file_name_base <- "Stacked barplot"
# 
# # Save to PDF
# ggsave(
#   filename = paste0(output_folder, file_name_base, ".pdf"),
#   plot = stacked_barplot_v2,
#   width = PUB_WIDTH * 2,
#   height = PUB_HEIGHT * 2,
#   units = "in",
#   device = cairo_pdf,  # ใช้ cairo_pdf เพื่อการแสดงผลที่ดีกว่าและรองรับอักขระพิเศษ
#   fallback_resolution = PUB_DPI
# )

# Microbiota plot
bracken_tibble = as.data.frame(bracken_matrix) %>% 
  rownames_to_column(var = "Sample ID") %>% 
  as_tibble()

bracken_tibble_location <- as.data.frame(bracken_tibble) %>% 
  left_join(
    sample_info %>% select(`Sample ID`, Location),
    by = "Sample ID")

# Check all tibble structure
glimpse(bracken_tibble_location)

# Data for ggbartax
data_for_ggbartax <- bracken_tibble_location %>% 
  arrange(Location) %>% 
  column_to_rownames(var = "Sample ID")

# Plot
pclass <- ggbartax(obj=data_for_ggbartax, facetNames= "Location",topn = 10) +
  facet_grid(~ Location, scales = "free_x", space = "free_x") +
  
  scale_y_continuous(expand = c(0, 0), limits = c(0, 101)) +  # ทำให้แท่งกราฟเริ่มที่ 0 พอดี และเผื่อที่ด้านบนเล็กน้อย
  scale_fill_manual(values=c(colorRampPalette(RColorBrewer::brewer.pal(12,"Set3"))(31))) +
  
  guides(fill= guide_legend(keywidth = 1.5, keyheight = 1.5, ncol=2))+
  
  labs(
    title = "Relative Abundance of Dominant Taxa",
    subtitle = "Top 14 most abundant genera are shown; others are grouped.",
    x = "Sample ID",
    y = "Relative Abundance"
  ) +
  
  theme_bw(base_size = 18) + # เพิ่ม base_size เพื่อให้ตัวอักษรโดยรวมใหญ่ขึ้น
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5, face = "italic"),
    
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), # ปรับมุม
    axis.title = element_text(face = "bold"),
    axis.line = element_line(linewidth = 0.4),
    axis.ticks = element_line(linewidth = 0.4),
    
    legend.title = element_text(face = "bold"),
    legend.background = element_blank(), # เอาพื้นหลัง legend ออก
    legend.key.size = unit(0.6, "cm"),
    
    panel.grid.major = element_blank(), # เอาเส้นกริดแนวตั้งออกเพื่อให้แท่งบาร์เด่นขึ้น
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", linewidth = 1),
    panel.spacing.x = unit(0.5, "lines"),   # ลดช่องว่างของ Facet
    
    # ปรับแต่งแถบหัวเรื่อง Facet
    strip.background = element_rect(fill = "gray85", color = "black"),
    strip.text.x = element_text(face = "bold", size = 14),
    
    text = element_text(family = "times_new_roman")
  )

print(pclass)

# Export for publication
output_folder <- "C:/Users/patch/Desktop/"
file_name_base <- "Stacked barplot"

# Save to PDF
ggsave(
  filename = paste0(output_folder, file_name_base, ".pdf"),
  plot = pclass,
  width = PUB_WIDTH * 2,
  height = PUB_HEIGHT * 2,
  units = "in",
  device = cairo_pdf,  # ใช้ cairo_pdf เพื่อการแสดงผลที่ดีกว่าและรองรับอักขระพิเศษ
  fallback_resolution = PUB_DPI
)


# Step 4: Stacked barplot with BWA and Kraken comparison ------------------



# Step 6: Boxplot ---------------------------------------------------------

# Kruskal test for significant difference test
res.kruskal_test <- alpha_diversity_df %>% kruskal_test(shannon_index ~ Location)
res.kruskal_test  # (if p > 0.05 means not significant)
# kruskal_statistic <- round(res.kruskal_test$statistic, 2) # ค่า Chi-squared

# If significant then perform wilcox test (post-hoc)
wilcox_test <- alpha_diversity_df %>% 
  wilcox_test(shannon_index ~ Location, p.adjust.method = "bonferroni") %>% 
  add_xy_position(x = "Location")
wilcox_test

# Extract only p-value adjusted
p_values_matrix <- wilcox_test %>%
  select(group1, group2, p.adj) %>%
  as.data.frame()

# Define p values name vector
p_values_named_vector <- p_values_matrix$p.adj
names(p_values_named_vector) <- paste(p_values_matrix$group1, p_values_matrix$group2, sep = "-")

# Convert p-values to full matrix
full_matrix <- multcompLetters(p_values_named_vector, compare= "<",
                               threshold = 0.05, Letters = letters)

# Create data frame with label "a, b, c" in each Location
df_label <- data.frame(
  Location = names(full_matrix$Letters),
  Label = full_matrix$Letters
)

# Find height for label locations
label_positions <- alpha_diversity_df %>%
  group_by(Location) %>%
  summarise(y_position = max(shannon_index)) %>%  # slightly higher
  left_join(df_label, by = "Location")

# พล็อตกราฟ (significant difference)
alpha_boxplot_v2 <- 
  ggplot(alpha_diversity_df, aes(x = Location, y = shannon_index, fill = Location)) +
  geom_violin(alpha = 0.7, color = NA) + # Violin plot เป็น optional ถ้าอยากให้กราฟดูคลีนขึ้น อาจจะเอาออก
  geom_boxplot(width = 0.6, alpha = 0.85, outlier.shape = NA) +  # ซ่อน outlier ของ boxplot เพื่อใช้ jitter แทน
  geom_jitter(
    aes(color = Location, fill = Location),
    width = 0.15,
    size = 2.5,
    stroke = 0.5,   # เพิ่มความหนาเส้นขอบ
    # shape = 21,   # ใช้ shape 21 ที่สามารถใส่สีพื้น (fill) และสีขอบ (color) ได้
    show.legend = FALSE
  ) +
  geom_text(
    data = label_positions,
    aes(x = Location, y = y_position, label = Label),
    size = 5,
    # fontface = "bold",
    vjust = -1,
    hjust = -1,
    family = "times_new_roman"
  ) +
  
  # stat_pvalue_manual(
  #   wilcox_test,
  #   label = "p.adj.signif", # แสดงเป็นสัญลักษณ์ (ns, *, **, ***)
  #   tip.length = 0.01,
  #   size = 8,
  #   bracket.size = 0.8,
  #   family = "times_new_roman"
  # ) +
  
  # scale_fill_manual(values = location_palettes$light) +
  scale_colour_manual(values = c("#c56761", "#968512", "#00952d",
                                 "#129c9f", "#4e7dcc", "#c450b6")
  ) +
  
  # ตรึงแกน y ที่ 0,0
  scale_y_continuous(limits = c(0, max(label_positions$y_position) + 0.4),
                     expand = c(0, 0)
  ) +
  
  # guides(
  #   fill = guide_legend(title = "Location",
  #                       override.aes = list(size = 8, colour = "black")),
  #   color = "none" # ซ่อน Legend ของสีขอบ (เพราะมันซ้ำกับสีพื้น)
  # ) +
  
  labs(
    title = "Alpha Diversity Comparison by Location",
    subtitle = paste("Kruskal-Wallis, p =", format.pval(res.kruskal_test$p, digits = 3)),
    x = "Location",
    y = "Shannon Diversity Index"
  ) +
  
  theme_bw(base_size = 16) + # เพิ่ม base_size เพื่อให้ตัวอักษรโดยรวมใหญ่ขึ้น
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5, size = 16, face = "italic"),
    
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    axis.title = element_text(face = "bold"),
    axis.line = element_line(size = 0.4),
    axis.ticks = element_line(size = 0.4),
    
    legend.title = element_text(face = "bold"),
    legend.background = element_blank(), # เอาพื้นหลัง legend ออก
    
    panel.background = element_rect(fill = "#ebebeb"),
    panel.grid.major = element_line(color = "white"),
    panel.grid.minor = element_line(color = "white"),
    panel.border = element_rect(color = "black", linewidth = 1.2),
    
    
    text = element_text(family = "times_new_roman")
  )

print(alpha_boxplot_v2)

# Export for publication
output_folder <- "C:/Users/patch/Desktop/"
file_name_base <- "Boxplot"

# Save to PDF
ggsave(
  filename = paste0(output_folder, file_name_base, ".pdf"),
  plot = alpha_boxplot_v2,
  width = PUB_WIDTH * 1,
  height = PUB_HEIGHT * 1.25,
  units = "in",
  device = cairo_pdf,  # ใช้ cairo_pdf เพื่อการแสดงผลที่ดีกว่าและรองรับอักขระพิเศษ
  fallback_resolution = PUB_DPI
)

# Non-significant difference case
alpha_boxplot_nonsig <- 
  ggplot(alpha_diversity_df, aes(x = Location, y = shannon_index, fill = Location))+
  geom_violin(alpha = 0.7, color = NA) + # Violin plot เป็น optional ถ้าอยากให้กราฟดูคลีนขึ้น อาจจะเอาออก
  geom_boxplot(width = 0.6, alpha = 0.85, outlier.shape = NA) +
  geom_jitter(
    aes(color = Location), # ไม่ต้องมี fill = Location ซ้ำซ้อน
    width = 0.15,
    size = 2.5,
    stroke = 0.5,
    show.legend = FALSE
  ) +
  
  # scale_fill_manual(values = location_palettes$light) +
  scale_colour_manual(values = c("#c56761", "#968512", "#00952d",
                                 "#129c9f", "#4e7dcc", "#c450b6")
  ) +
  
  scale_y_continuous(
    limits = c(0, NA), 
    expand = expansion(mult = c(0, 0.1)) # เพิ่มพื้นที่ด้านบน 10%
  ) +
  
  guides(
    fill = guide_legend(title = "Location"),
    color = "none"
  ) +
  
  labs(
    title = "Alpha Diversity Comparison by Location",
    subtitle = paste("Kruskal-Wallis test, p =", format.pval(res.kruskal_test$p, digits = 3)),
    x = "Location",
    y = "Shannon Diversity Index"
  ) +
  
  theme_bw(base_size = 16) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5, size = 16, face = "italic"),
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    axis.title = element_text(face = "bold"),
    axis.line = element_line(size = 0.4),
    axis.ticks = element_line(size = 0.4),
    legend.title = element_text(face = "bold"),
    legend.background = element_blank(),
    panel.background = element_rect(fill = "#ebebeb"),
    panel.grid.major = element_line(color = "white"),
    panel.grid.minor = element_line(color = "white"),
    panel.border = element_rect(color = "black", linewidth = 1.2),
    text = element_text(family = "times_new_roman")
  )

print(alpha_boxplot_nonsig)

# Export for publication
output_folder <- "C:/Users/patch/Desktop/"
file_name_base <- "Boxplot_NonSignificant"

# Save to PDF
ggsave(
  filename = paste0(output_folder, file_name_base, ".pdf"),
  plot = alpha_boxplot_nonsig,
  width = PUB_WIDTH * 1,
  height = PUB_HEIGHT * 1.25,
  units = "in",
  device = cairo_pdf,
  fallback_resolution = PUB_DPI
)


# Step 7: PCoA plot --------------------------------------------------------

# Beta Diversity Analysis & PCoA
beta_diversity_dist <- vegdist(bracken_matrix, method = "bray") # ใช้ bracken matrix จาก  Step 4: Boxplot ได้เลย
pcoa_result <- wcmdscale(beta_diversity_dist, k = 2, eig = TRUE)

# ดึงค่า variance explained ออกมา
variance_explained <- pcoa_result$eig / sum(pcoa_result$eig)

# สร้างข้อความ label ให้กับแกนแต่ละแกน
pcoa_axis1_lab <- sprintf("PCoA 1 (%.2f%%)", variance_explained[1] * 100)
pcoa_axis2_lab <- sprintf("PCoA 2 (%.2f%%)", variance_explained[2] * 100)

# สร้าง dataframe จาก pcoa_result
pcoa_df <- as.data.frame(pcoa_result$points) %>% 
  rename("PCoA1" = Dim1, "PCoA2" = Dim2) %>%   # เปลี่ยนชื่อ column
  rownames_to_column(var = "sample_id") %>%    # เปลี่ยน rownames ให้เป็น column ชื่อ sample_id
  left_join(metadata_df, by = "sample_id")  # join ข้อมูล metadata จาก alpha_diversity_df

# ทดสอบนัยสำคัญทางสถิติด้วย PERMANOVA
set.seed(1024)
permanova_result <- adonis2(beta_diversity_dist ~ Location, data = pcoa_df)

# สร้าง label เพื่อแสดงผลบนกราฟ 
permanova_label <- sprintf("PERMANOVA\nR² = %.3f, p = %.3f",
                           permanova_result$R2[1],
                           permanova_result$`Pr(>F)`[1])

# Check outlier ที่ Check_Outlier.R
# สร้าง Label ให้ Outlier
outlier_id <- pcoa_df$sample_id[7]

# สร้างคอลัมน์ "label" ขึ้นมาใหม่
pcoa_df <- pcoa_df %>% 
  mutate(label = ifelse(sample_id == outlier_id, sample_id, ""))

# พล็อตกราฟ
pcoa_plot_v2 <-
  ggplot(pcoa_df, aes(x = PCoA1, y = PCoA2, fill = Location,
                      shape = Location, text = paste0("Sample ID:", sample_id))
  ) +
  
  # เพิ่มเส้นประที่จุด 0
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray70") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray70") +
  
  # เพิ่ม Confidence Elipse เพื่อแสดงการจัดกลุ่ม
  # stat_ellipse(aes(group = Location), level = 0.95, geom = "polygon", alpha = 0.15, show.legend = FALSE) +
  
  # เพิ่มจุดของข้อมูล
  geom_point(shape = 21,
             size = 4,
             stroke = 1.2,
             alpha = 0.8,
             color = "black"
             ) +
  
  # เพิ่ม label ในจุดที่ drift ออกมา
  # geom_text_repel(
  #   aes(label = label),
  #   color = "black", size = 4,
  #   box.padding = 0.5, point.padding = 0.5,
  #   max.overlaps = Inf,
  #   family = "times_new_roman"
  # ) +
  
  # เพิ่ม label ผลสถิติที่มุมขวาล่าง
  # annotate("text", x = Inf, y = -Inf, label = permanova_label,
  #          hjust = 1.1, vjust = -0.5, size = 5, family = "times_new_roman",
  #          fontface = "plain") +
  
  stat_chull(
    aes(fill = Location, color = Location), # บอกให้สีของ "เส้น" เปลี่ยนตาม Location
    alpha = 0.2,
    geom = "polygon",
    show.legend = FALSE,    # ไม่ต้องแสดง legend ของเส้นซ้ำซ้อน
  ) +
  
  scale_fill_manual(values = c("#c56761", "#968512", "#00952d",
                               "#129c9f", "#4e7dcc", "#c450b6"), 
                    name = "Location") +
  scale_color_manual(values = c("#c56761", "#968512", "#00952d",
                                "#129c9f", "#4e7dcc", "#c450b6"),
                     name = "Location") + # ใช้สีชุดเดียวกันสำหรับเส้น ellipse
  
  guides(
    fill = guide_legend(title = "Location", override.aes = list(shape = 21, size = 8, linetype = 0, color = "black")),
    color = "none" # ซ่อน Legend ของสีขอบ (เพราะมันซ้ำกับสีพื้น)
  ) +
  
  labs(
    title = "Beta Diversity based on Bray-Curtis Dissimilarity",
    subtitle = "Principal Coordinate Analysis (PCoA)",
    x = pcoa_axis1_lab,
    y = pcoa_axis2_lab
  ) +
  
  theme_bw(base_size = 18) + # เพิ่ม base_size เพื่อให้ตัวอักษรโดยรวมใหญ่ขึ้น
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    axis.title = element_text(face = "bold"),
    axis.line = element_line(linewidth = 0.4),
    axis.ticks = element_line(linewidth = 0.4),
    
    legend.title = element_text(face = "bold"),
    legend.background = element_blank(), # เอาพื้นหลัง legend ออก
    
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", linewidth = 1.2),
    
    text = element_text(family = "times_new_roman")
  )

print(pcoa_plot_v2)

# Export for publication
output_folder <- "C:/Users/patch/Desktop/"
file_name_base <- "PCoA_plot"

# Save to PDF
ggsave(
  filename = paste0(output_folder, file_name_base, ".pdf"),
  plot = pcoa_plot_v2,
  width = PUB_WIDTH * 1.5,
  height = PUB_HEIGHT * 1.5,
  units = "in",
  device = cairo_pdf,  # ใช้ cairo_pdf เพื่อการแสดงผลที่ดีกว่าและรองรับอักขระพิเศษ
  fallback_resolution = PUB_DPI
)


# Step 8: Rarefaction Curve -------------------------------------------------------

# เตรียม OTU Table ให้ถูกต้อง
otu_tab_for_mpse <- bracken_data_wide %>% 
  column_to_rownames(var = "sample_id") %>% 
  t() %>% 
  as.data.frame()
# setNames(str_remove(colnames(.), "_\\d+$"))  # ให้ colnames ของ . (pipe ณ ปัจจุบัน) ตรงกับ rownames ของ sample_info_for_mpse

# เตรียม Sample Data ให้ถูกต้อง
sample_info_for_mpse <- alpha_diversity_df %>% 
  select(sample_id, Location) %>% 
  distinct() %>% 
  column_to_rownames(var = "sample_id") %>%  # ทำให้ column ที่ชื่อ join_id เป็น rownames
  .[colnames(otu_tab_for_mpse), , drop = FALSE] # จัดลำดับ rownames ของ Sample Data ให้ตรงกับ colnames ของ OTU Table

# สร้าง MPSE Object
mpse_obj <- MPSE(
  assays = list(Abundance = otu_tab_for_mpse),
  colData = sample_info_for_mpse
)

# คำนวณข้อมูลสำหรับ Rarefaction Curve
set.seed(1024)
mpse_rare <- mp_cal_rarecurve(
  .data = mpse_obj,
  .abundance = Abundance,
  .group = "Location", # บอกให้จัดกลุ่มตาม Location!
  chunks = 400 # เพิ่มความละเอียดของเส้นโค้ง
)

# ดึง list ของ tibble ออกมาจาก object
rare_data_list <- mpse_rare@colData$RareAbundanceRarecurve  # มี 11 list ใหญ่

# รวม list ของ tibble  ทั้งหมดให้กลายเป็น dataframe ใหญ่
rare_data_tidy <- bind_rows(rare_data_list, .id = "Sample") # เพิ่มคอลัมน์ Sample เข้าไปด้วย

# พล็อตกราฟเฉพาะดัชนี Observe
rarefaction_plot_v2 <- 
  ggplot(
    
    # กรองข้อมูลให้เหลือเฉพาะ Alpha == "Observe"
    filter(rare_data_tidy, Alpha == "Observe"), aes(x = readsNums, y = value, group = Sample, color = Location)) +
  geom_line(linewidth = 0.8, alpha = 0.8) +
  
  scale_color_manual(values = location_palettes$light) +
  
  labs(
    title = "Rarefaction Curve by Location",
    x = "Sequencing Depth (Reads)",
    y = "Observed Features (Richness)",
    color = "Location"
  ) +
  
  theme_bw(base_size = 16) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.title = element_text(face = "bold"),
    legend.title = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    text = element_text(family = "times_new_roman")
  )

print(rarefaction_plot_v2)

# พล็อตหลาย ๆ ดัชนี
rarefaction_plot_multiple_v2 <- 
  ggplot(rare_data_tidy, aes(x = readsNums, y = value, group = Sample, color = Location)) +
  geom_line(linewidth = 0.7, alpha = 0.7) +
  
  facet_wrap(~ Alpha, scales = "free_y", ncol = 3) +  # ให้แกน Y ของแต่ละ Panel มีสเกลของตัวเอง
  
  # scale_color_manual(values = location_palettes$light) +
  
  labs(
    title = "Rarefaction Curves for Multiple Diversity Indices",
    x = "Sequencing Depth (Reads)",
    y = "Alpha Diversity Value", # ใช้ชื่อกลางๆ เพราะแต่ละ Panel มีดัชนีของตัวเอง
    color = "Location"
  ) +
  
  theme_bw(base_size = 16) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.title = element_text(face = "bold"),
    legend.title = element_text(face = "bold"),
    legend.background = element_blank(),
    panel.grid.major = element_line(color = "grey90"), # อาจจะเปิด Grid อ่อนๆ ไว้เพื่อให้ดูง่าย
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", linewidth = 1),
    strip.background = element_rect(fill = "grey85", color = "black"), # ปรับแต่งหัวข้อของแต่ละ Panel
    strip.text = element_text(face = "bold"),
    text = element_text(family = "times_new_roman")
  )

print(rarefaction_plot_multiple_v2)

# คำนวณเส้นค่าเฉลี่ยแต่ละกลุ่ม
rare_data_summary <- rare_data_tidy %>% 
  group_by(Location, Alpha, readsNums) %>% 
  summarise(
    mean_value = mean(value),
    se_value = sd(value) / sqrt(n()),
    .groups = "drop"
  ) %>% 
  filter(!is.na(se_value))

# เตรียมข้อมูลสำหรับ label (เลือกเฉพาะจุดสุดท้ายของแต่ละเส้น)
label_data <- rare_data_summary %>%
  filter(Alpha %in% c("Observe", "Shannon")) %>%
  group_by(Location) %>%
  filter(readsNums == max(readsNums))

# Hybrid plot
rarefaction_plot_final <- 
  ggplot(filter(rare_data_summary, Alpha %in% c("Observe", "Shannon")),
         aes(x = readsNums, y = mean_value, color = Location, fill = Location)) +
  geom_ribbon(aes(ymin = pmax(mean_value - se_value), ymax = mean_value + se_value),
              alpha = 0.15, color = NA) +
  geom_line(linewidth = 1.2) +  # เส้นค่าเฉลี่ย
  # geom_smooth(method = "loess", se = FALSE, linewidth = 1.3, span = 0.4) +  # <<<<< เส้นเรียบ
  
  geom_text_repel(
    data = label_data,
    aes(label = Location),
    color = "black",
    nudge_x = 50, # ขยับ label ออกไปทางขวา
    direction = "both",
    hjust = 0,
    segment.size = NA,
    size = 4.5,
    box.padding = 0.4  
  ) +
  
  facet_wrap(~ Alpha, scales = "free_y", ncol = 2) +
  
  # scale_color_manual(values = location_palettes$light) +
  scale_colour_manual(values = c("#c56761", "#968512", "#00952d",
                                 "#129c9f", "#4e7dcc", "#c450b6")
  ) +
  
  guides(color = guide_legend(nrow = 2, byrow = TRUE)) +
  scale_x_continuous(
    breaks = seq(0, 2500, by = 500), # กำหนดให้มีขีดที่ 0, 500, ..., 2500
    limits = c(0, 2850),             # ขยายขอบเขตการวาดจริงไปถึง 2850 (ปรับเลขได้)
    expand = c(0, 0)                 # ทำให้กราฟเริ่มที่ 0 พอดี ไม่มีช่องว่าง
  ) +
  
  labs(
    title = "Rarefaction Curve by Location",
    x = "Sequencing Depth (Reads)",
    y = "Alpha Diversity Value",
    color = "Location"
  ) +
  
  theme_bw(base_size = 16) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.title = element_text(face = "bold"),
    legend.title = element_text(face = "bold"),
    legend.background = element_blank(),
    legend.position = "none",
    # panel.grid.major = element_line(color = "grey90"),
    # panel.grid.minor = element_blank(),
    # panel.border = element_rect(color = "black", linewidth = 1),
    
    panel.background = element_rect(fill = "#ebebeb"),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.border = element_rect(color = "black", linewidth = 1.2),
    
    strip.background = element_rect(fill = "grey85", color = "black"), # <<<< ปรับแต่งหัวข้อ Panel
    strip.text = element_text(face = "bold", color = "black"), # <<<< ปรับแต่งตัวอักษรหัวข้อ Panel
    
    text = element_text(family = "times_new_roman")
  )

print(rarefaction_plot_final)

# Export for publication
output_folder <- "C:/Users/patch/Desktop/"
file_name_base <- "Rarefaction_Curve"

# Save to PDF
ggsave(
  filename = paste0(output_folder, file_name_base, ".pdf"),
  plot = rarefaction_plot_final,
  width = PUB_WIDTH * 2.75,
  height = PUB_HEIGHT * 1.5,
  units = "in",
  device = cairo_pdf,  # ใช้ cairo_pdf เพื่อการแสดงผลที่ดีกว่าและรองรับอักขระพิเศษ
  fallback_resolution = PUB_DPI
)


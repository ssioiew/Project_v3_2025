# ------------------------------------------------------
# File: summarize_bracken_data.R
# Purpose: To consolidate and summarize microbial abundance data from Bracken reports.
# Author: Pirada Naewkam
# Date: 2025-08-27
#
# Description:
#    This script reads all Bracken report files from specified directories
#    (based on different trimming conditions and analysis methods). It then
#    merges the data into a single long-format table, pivots it into a wide format,
#    and filters for the top 15 most abundant organisms before saving the
#    final table to an Excel file.
#
# Input: Bracken report files (.bracken.G.report) from different folders.
# Output: An Excel file named "Test_R_Code_Summary_Table.xlsx" containing the summarized data.
# ------------------------------------------------------


# Step 0: Setup All Environments ------------------------------------------

# Load All Libraries
library(tidyverse)
library(writexl)

# Setup All Paths
base_path <- "Kraken Data"
data_sources <- tribble(
  ~condition, ~method, ~folder_path,
  "not trim", "Kraken", file.path(base_path, "not trim/Kraken"),
  "not trim", "Map", file.path(base_path, "not trim/Map"),
  "fastp", "Kraken", file.path(base_path, "fastp/Kraken"),
  "fastp", "Map", file.path(base_path, "fastp/Map"),
  "trimmomatic", "Kraken", file.path(base_path, "trimmomatic/Kraken"),
  "trimmomatic", "Map", file.path(base_path, "trimmomatic/Map")
)


# Step 1: Data Wraggling --------------------------------------------------

# สร้างฟังก์ชันสำหรับอ่านและประมวลผลไฟล์
read_bracken_file <- function(filepath, condition_name, method_name) {
  
  # ดึงชื่อ Sample_ID จากชื่อไฟล์
  sample_id <- basename(filepath) %>% 
    str_remove("\\.bracken\\.G\\.report$")  # เอานามสกุลของไฟล์ออก (เหลือแต่ชื่อ Sample_ID)
  
  # อ่านไฟล์ Bracken แบบ tab-seperated
  bracken_data <- read_tsv(filepath, show_col_types = FALSE) %>% 
    select(organism = name, fraction = fraction_total_reads) %>% 
    mutate(
      percentage = fraction * 100, # แปลง fraction เป็น %
      sample = sample_id,
      condition = condition_name,
      method = method_name
    ) %>% 
    select(sample, condition, method, organism, percentage) # จัดลำดับ column
  
  return(bracken_data)
}

# วนลูปไฟล์ทั้งหมดรวมเป็นตารางเดียว ด้วย pmap (วนลูปตาม data_sources)
all_data_list <- list() # สร้าง list ว่างรอเก็บข้อมูล

all_data_list <- pmap(data_sources, function(condition, method, folder_path) {
  
  # หา bracken reports ทั้งหมดในโฟลเดอร์นั้น ๆ
  file_to_read <- list.files(
    path = folder_path,
    pattern = "\\.bracken\\.G\\.report$",
    full.names = TRUE
  )
  
  # ตรวจสอบจำนวนไฟล์
  if (length(file_to_read) == 0) {
    warning(paste("No files found in:", folder_path))
    return(NULL)
  }
  
  # วนอ่านทุกไฟล์ใน list แล้วรวมเป็น Dataframe เดียว
  map_df(file_to_read, ~read_bracken_file(., condition_name = condition, method_name = method))
})

# รวมข้อมูลจากทุก ๆ folder เข้าด้วยกันเป็นตารางใหญ่ตารางเดียว (long format)
long_format_data <- bind_rows(all_data_list)


# Step 2: Convert Data to Wide Format -------------------------------------

final_table <- long_format_data %>% 
  mutate(
    condition = factor(condition, levels = c("not_trim", "fastp", "trimmomatic")),
    method = factor(method, levels = c("Kraken", "Map"))
  ) %>% 
  arrange(sample, condition, method) %>%  # จัดเรียงข้อมูลก่อน pivot
  pivot_wider(
    names_from = c(condition, method, sample),   # สร้าง column จากสองตัวนี้
    values_from = percentage,
    values_fill = 0.0 # ถ้าไม่มีข้อมูลให้ใส่ 0
  ) %>%
  arrange(organism) # จัดเรียงแถวตามชื่อ organism


# Step 3: Filter for Top 15 Organisms ---------------------------------------

# ใช้ threshold
# กรองเอาเฉพาะ % ของ organism ที่มีค่าสูงใน sample ใด sample หนึ่ง
# top_organisms_table <- long_format_data %>% 
#   group_by(organism) %>% 
#   summarise(max_percentage = max(percentage)) %>% 
#   filter(max_percentage > 50) %>% 
#   pull(organism)

# นำ list top organisms ไปกรอกตารางสุดท้าย
# final_table_top_only <- final_table %>% 
#   filter(organism %in% top_organisms_table)

# บันทึกผล
# write_xlsx(final_table_top_only, path = "Test_R_Code_Summary_Table.xlsx")

# ใช้ Top N
# organism_max_abundance <- long_format_data %>% 
#   filter(!str_detect(organism, "Homo")) %>% 
#   group_by(organism) %>% 
#   summarise(max_percentage = max(percentage), .groups = "drop")
# 
# top15_organisms_list <- organism_max_abundance %>% 
#   arrange(desc(max_percentage)) %>% 
#   slice_head(n = 15) %>% 
#   pull(organism)

# เพิ่ม host เข้าไปด้วย
# top15_with_host <- c("Homo", top15_organisms_list)

# นำ list top15 ไปกรองในตารางสุดท้าย
# final_table_top15 <- final_table %>% 
#   filter(organism %in% top15_with_host)

# บันทึกผล
# write_xlsx(
#   list(   # แยก sheet ใน excel
#     Top_15_Organisms = final_table_top15,
#     All_Organisms = final_table
#   ),
#   path = "Organism_Summary_Top15.xlsx"
# )

summary_composition <- long_format_data %>% 
  mutate(
    category = if_else(str_detect(organism, "Homo"), "Host (Homo)", "Microbe")
  ) %>% 
  group_by(sample, condition, method, category) %>% 
  summarise(total_percentage = sum(percentage), .groups = "drop")

# เปลี่ยนเป็น wide format และคำนวณ unclassified เผื่อไว้ด้วย
composition_table_wide <- summary_composition %>% 
  pivot_wider(
    names_from = category,
    values_from = total_percentage,
    values_fill = 0.0 # ถ้าไม่มีให้ใส่ 0 (ทั้ง Host & Microbes)
  ) %>% 
  
  # คำนวณที่เหลือเป็น unclassified
  mutate(
    Unclassified = pmax(0, 100 - (`Host (Homo)` + Microbe))) %>% 
  select(sample, condition, method, `Host (Homo)`, Microbe, Unclassified) %>%
  arrange(condition, method, sample)

# บันทึกไฟล์
write_xlsx(
  composition_table_wide,
  path = "Sample_Composition_Summary02.xlsx"
)


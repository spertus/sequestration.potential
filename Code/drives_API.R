library(drivesR)
# Script to generate DRIVES API key

directus_token <- read_directus_pat("DRIVES_PAT.txt")

set_default_token(directus_token)

all_table_names <- options("drivesR.default.tablevec")[[1]]

# downloads all tables
db_tables <- import_db_tables(tablevec = all_table_names, savedir = "Data/DRIVES", savename = "all_drives_data")


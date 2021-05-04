#
# Autor(es): Camilo Avellaneda
# Mantenimiento: Camilo Avellaneda
# Fecha creación: 04/05/2021
#==============================================
getwd()
if(!("input" %in% list.files())){dir.create("input")}
if(!("output" %in% list.files())){dir.create("output")}
if(!("src" %in% list.files())){dir.create("src")}
if(!("docs" %in% list.files())){dir.create("docs")}
if(!("reports" %in% list.files())){dir.create("reports")}

require(pacman)
p_load(tidyverse,openxlsx,here,haven,data.table,fda,fda.usc,roahd,fdaoutlier,caret)
p_load(maptools,sgeostat,scatterplot3d,car,fields,gstat,geoR)

options(scipen=999)
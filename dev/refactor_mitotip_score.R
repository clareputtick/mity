mtip <- read.csv("mitylib/annot/mitotip_score_fixed_del.csv")
mtip$MitoTip_interpretation <- "likely benign"
mtip$MitoTip_interpretation[orig$MitoTip_score > 8.44] <- "possibly benign"
mtip$MitoTip_interpretation[orig$MitoTip_score > 12.66] <- "possibly pathogenic"
mtip$MitoTip_interpretation[orig$MitoTip_score > 16.25] <- "likely pathogenic"

# write.csv(mtip, "dev/mitotip_score_fixed_del_fixed.csv", row.names=FALSE)
write.csv(mtip, "mitylib/annot/mitotip_score_fixed_del.csv", row.names=FALSE)
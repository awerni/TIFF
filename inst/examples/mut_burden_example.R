
sql1 <- paste0("select count(*), tissuename from tissue.processedsequenceview WHERE tissuename in (select tissuename FROM tissue.processedsequenceextended ",
               "WHERE tissuename in (select tissuename from tissue.tissue where tumortype = 'stomach adenocarcinoma') and enst in ",
               "(select enst from gene natural join transcript where symbol  in ('POLE', 'POLD1') and iscanonical) AND coarse(aamutation) = 'wt') group by tissuename")

data_wt <- getPostgresql(sql1)


sql2 <- paste0("select count(*), tissuename from tissue.processedsequenceview WHERE tissuename in (select tissuename FROM tissue.processedsequenceextended ",
               "WHERE tissuename in (select tissuename from tissue.tissue where tumortype = 'stomach adenocarcinoma') and enst in ",
               "(select enst from gene natural join transcript where symbol  in ('POLE', 'POLD1') and iscanonical) AND coarse(aamutation) = 'mut') group by tissuename")

data_mut <- getPostgresql(sql2)

summary(data_mut$count)
summary(data_wt$count)
args = commandArgs(trailingOnly = TRUE)

p_value = read.delim(file=args[1],header = T,sep = '\t')
p_12 = as.double(p_value$p_12)
p_78 = as.double(p_value$p_78)

# FDR correction
p_12_corr = p.adjust(p_12,method='BH')
p_78_corr = p.adjust(p_78,method='BH')

# out = data.frame(column1=p_value$cluster,column2=p_12_corr,column3=p_78_corr)
# colnames(out) = c('cluster','p_12_adj','p_78_adj')
# write.table(out,file=args[2],sep='\t',row.names=F,quote = F)

# output
out = data.frame(column1=p_value$cluster,column2=p_12,column3=p_78,columns4 = p_12_corr,columns5 = p_78_corr)
colnames(out) = c('cluster','p_12','p_78','p_12_adj','p_78_adj')
write.table(out,file=args[2],sep='\t',row.names=F,quote = F)
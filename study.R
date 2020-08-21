plot1 <- plot_allyears_allcompgrid()
ggsave(paste0('allcomponents_', country_code,'.png'), plot = plot1, width = 6, height = 6)

p2 <- plot_Ex_by_comp()
ggsave(paste0('comp_imp_', country_code,'.png'), plot = p2, width = 8, height = 8)


#QX plots
plot2 <- plot_allyears_qx()
ggsave(paste0('allyearsQx_', country_code,'.png'), plot = plot2, width = 6, height = 6)
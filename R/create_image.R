# Plot image of sparse matrix

create_image = function(matrix){
  
  matrix %>% as.vector %>% 
    tibble(value = ., row = rep(1:nrow(matrix), times = ncol(matrix)),
           col = rep(1: ncol(matrix), each = nrow(matrix))) %>%
    ggplot(aes(x = col, y = row, shape = value)) + scale_x_continuous(position="top") +
    geom_point(size = 2) + ylab("row") + scale_y_reverse() + 
    scale_shape_manual(name = "", values = c(16,1)) +
    theme_minimal() # + annotate(geom="text", x=48, y=94, label = expression("Group lasso (" * omega[n] * " = 0.670)"))
  
}

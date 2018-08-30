#Cutoff graph

Cells_included <- NULL; Cutoff_percent <- seq(from = 0, to = 1, by = 0.01)

for (i in Cutoff_percent){
  Cell_number <- ConPoly[(ConPoly$Area >= i),]
  Cell_number <- length(Cell_number)
  Cells_included <- append(Cells_included, Cell_number)}
Attrition <- data.frame(Cutoff_percent, Cells_included)
Attrition$Percent_compared_to_Sarah <- 100*Cells_included/max(Cells_included)

#Percent
ggplot(Attrition, aes(x=Cutoff_percent*100, y=Percent_compared_to_Sarah))+geom_point()+
  #labs(x = "Cutoff percent", y = "Percent of cells (compared to original)")+
  scale_y_continuous(name = "Percent of cells remaining", breaks = c(seq(from = 35, to = 100, by = 10),100))+
  scale_x_continuous(name = "Acceptable new cell size (percent compared to pre-clip)", breaks = seq(from = 0, to = 100, by = 10))+
  
  theme(
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 12),
    axis.line = element_line("black"),
    legend.position="none",
    #legend.title = element_blank()
    panel.background = element_blank()
  )

ggsave("./BC_ConnectivityProject/Cell cutoff graph - Percent.png", width = 10, height = 6)



#Total cells
ggplot(Attrition, aes(x=Cutoff_percent*100, y=Cells_included))+geom_point()+
  #labs(x = "Cutoff percent", y = "Percent of cells (compared to original)")+
  scale_y_continuous(name = "Number of cells remaining", breaks = seq(from = 350, to = 950, by = 100))+
  scale_x_continuous(name = "Acceptable new cell size (percent compared to pre-clip)", breaks = seq(from = 0, to = 100, by = 10))+
  
  theme(
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 12),
    axis.line = element_line("black"),
    legend.position="none",
    #legend.title = element_blank()
    panel.background = element_blank()
  )

ggsave("./BC_ConnectivityProject/Cell cutoff graph - Values.png", width = 10, height = 6)


rm(Cell_number, Cells_included, Cutoff_percent, Attrition)

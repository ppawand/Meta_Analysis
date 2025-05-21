# rank the predictor based on their relative predictive power

var_importance <- function (data) {
  cforest <- cforest(LRR ~ env.type + pH.level + ecosystem + irrigation.status +
                       magnitude.of.warming + warming.duration + Warming.technique,
                     data = data,
                     controls = (cforest_control(ntree = 1000, replace = T)))
  
  varimp <- varimp(cforest)
  plot(varimp)
  ranked.importance <- as.data.frame(sort(varimp, decreasing = T))
  names(ranked.importance) <- "index"
  ranked.importance$predictor <- row.names(ranked.importance)
  ranked.importance$predictor <- factor(ranked.importance$predictor,
                                        levels = ranked.importance$predictor[order(ranked.importance$index)])
  ranked_plot <- ggplot(ranked.importance, aes(x = predictor, y = index)) +
    geom_bar(stat = "identity", fill = "black", col = "black") +
    coord_flip() +
    theme_bw()
  return(ranked_plot)
}
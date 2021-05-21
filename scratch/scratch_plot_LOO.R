ploo <- gp$pred_LOO(se.fit = T)
gp$Z
loodf <- cbind(ploo, Z=gp$Z)
loodf
loodf$upper <- loodf$fit + 1.96 * loodf$se.fit
loodf$lower <- loodf$fit - 1.96 * loodf$se.fit
ggplot(loodf, aes(fit, Z)) +
  stat_smooth() +
  geom_abline(slope=1, intercept=0, color="red") +
  geom_segment(aes(x=lower, xend=upper, yend=Z), color="green") +
  geom_point()
# Add text with coverage, R-sq
coveragevec <- with(loodf, upper >= Z & lower <= Z)
coverage <- mean(coveragevec)
coverage
rsq <- with(loodf, 1 - (sum((fit-Z)^2)) / (sum((mean(Z)-Z)^2)))
rsq
ggplot(loodf, aes(fit, Z)) +
  stat_smooth() +
  geom_abline(slope=1, intercept=0, color="red") +
  geom_segment(aes(x=lower, xend=upper, yend=Z), color="green") +
  geom_point() +
  # geom_text(x=min(loodf$fit), y=max(loodf$Z), label="abc") +
  geom_text(x=-Inf, y=Inf, label=paste("Coverage:", signif(coverage,5)), hjust=0, vjust=1) +
  geom_text(x=-Inf, y=Inf, label=paste("R-sq:        ", signif(rsq,5)), hjust=0, vjust=2.2) +
  # geom_text(x=Inf, y=-Inf, label="def", hjust=1, vjust=0)
  xlab("Predicted values (fit)") +
  ylab("Actual values (Z)") +
  ggtitle("Calibration of leave-one-out (LOO) predictions")

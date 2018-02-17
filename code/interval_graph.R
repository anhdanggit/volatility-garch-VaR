library(ggplot2)
library(gridExtra)

mu.interval = read.csv("results/interval_mu.csv")
mu.interval$n = factor(mu.interval$n)
mu.interval$model = c("100.std", "1000.std", "100.norm", "1000.norm")
mu.interval$est_lower = mu.interval$est - mu.interval$sd
mu.interval$est_upper = mu.interval$est + mu.interval$sd

p = ggplot(mu.interval, aes(x=model, y=est, ymin=est-sd, ymax=est+sd)) + 
  geom_linerange(size=8, colour="#a6d8f0") +
  geom_hline(aes(x=0, yintercept=2), lty=1) + # change hline value
  geom_point(size=3, shape=21, fill="#008fd5", colour = "white", stroke = 1) +
  coord_flip() +
  ggtitle("theta - design 2") + # change tital
  theme_minimal()  

p

p5 = p

grid.arrange(p1, p2, p3, p4, p5, ncol=5)


# theoretical densities of norm and std
xseq<-seq(-4,4,.01)
densities<-dnorm(xseq, 0,1)
densities2 <- dt(xseq, df=10)
df = data.frame(xseq = xseq, den = densities, den2 = densities2)

ggplot(df, aes(x=xseq, y=densities)) +
  geom_line() +
  geom_line(aes(y=densities2), color="red")


tab <- read.csv("results/eq1_norm.csv")
stargazer::stargazer(tab)

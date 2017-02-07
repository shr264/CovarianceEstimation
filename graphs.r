library(gRbase)
library(Rgraphviz)
ug0 <- ug(~1:2+2:3+3:1+2:4+1:5)
ug0
pdf("graph2.pdf",width=3,height=3)
plot(ug0, attrs=list(node=list(label="foo", fillcolor="lightgreen")))
dev.off()
ug1 <- ug(~1:2+1:3+1:4+2:5+3:5)
pdf("graph3.pdf",width=3,height=3)
plot(ug1, attrs=list(node=list(label="foo", fillcolor="lightgreen")))
dev.off()
ag <- as(ug1, "matrix")
ag

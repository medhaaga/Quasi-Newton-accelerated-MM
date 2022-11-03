set.seed(1)
x <- seq(-2*pi, 2*pi, length=100)
pts1 <- runif(1)*(2*pi)
pts2 <- runif(1)*(2*pi)
pts3 <- pts2 + sin(pts2)
pts4 <- pts2 - sin(pts2)*((pts2-pts1)/(sin(pts2)-sin(pts1)))
pts5 <- pts2 - sin(pts2)*((pts3-pts2)/(sin(pts3)-sin(pts2)))

pdf(file = "secants.pdf", height=5, width = 7)
plot(x, sin(x), "l", xlab = "x", ylab = "G(x)")
points(x = c(pts1, pts2, pts3, pts4, pts5), y = c(sin(pts1), sin(pts2), sin(pts3), sin(pts4), sin(pts5)),  pch=19, cex=1.5)
abline(v = pi, lty=2)
segments(pts1, sin(pts1), pts2, sin(pts2), lty=3, col = "red", lwd=2)
segments(pts2, sin(pts2), pts3, sin(pts3), lty=2, col = "green3", lwd=2)
text(x = c(pts1-.3, pts2-.3, pts3-.3, pts4-.3, pts5-.3), y = c(sin(pts1), sin(pts2), sin(pts3), sin(pts4), sin(pts5)), c("A", "B", "C*", "C1", "C2"))
dev.off()

N <- 1000
start1 <- runif(N)^2*pi
start2 <- runif(N)^2*pi
itr1 <- rep(0, N)
itr2 <- rep(0, N)
conv1 <- rep(0, N)
conv2 <- rep(0, N)
tol <- 1e-7
for (i in 1:N){
  pts1 <- start1[i]
  pts2 <- start2[i]
  
  
  old <- pts1
  new <- pts2
  diff <- norm(new-old, "2")
  itr <- 0
  while(diff >= tol){
    foo <- new - sin(new)*((new-old)/(sin(new)-sin(old)))
    old <- new
    new <- foo
    diff <- norm(new-old, "2")
    itr <- itr+1
    }
  print(new)
  print(itr)
  itr1[i] <- itr
  conv1[i] <- new
  
  old <- pts2
  new <- pts2 + sin(pts2)
  diff <- 10
  tol <- 1e-7
  itr <- 0
  while(diff >= tol){
    foo <- old - sin(old)*((new-old)/(sin(new)-sin(old)))
    old <- foo
    new <- foo + sin(foo)
    diff <- norm(new-old, "2")
    itr <- itr+1
  }
  print(new)
  print(itr)
  itr2[i] <- itr
  conv2[i] <- new
  
}

print(quantile(itr1))
print(quantile(itr2))

pdf(file = "secants.pdf", height=5, width = 10)
par(mfrow = c(1,2))
pts1 <- pi/2 - 1.5
pts2 <- 1.5*pi - .1
pts3 <- pts2 + sin(pts2)
pts4 <- pts2 - sin(pts2)*((pts2-pts1)/(sin(pts2)-sin(pts1)))
pts5 <- pts2 - sin(pts2)*((pts3-pts2)/(sin(pts3)-sin(pts2)))

x <- seq(0, 2*pi, .01)
plot(x, sin(x), "l", xlab = "x", ylab = "G(x)", main = "Iteration 1")
points(x = c(pts1, pts2, pts3, pts4, pts5), y = c(sin(pts1), sin(pts2), sin(pts3), sin(pts4), sin(pts5)),  pch=19, cex=1.5)
abline(v = pi, lty=2)
segments(pts1, sin(pts1), pts2, sin(pts2), lty=3, col = "red", lwd=2)
segments(pts2, sin(pts2), pts3, sin(pts3), lty=2, col = "green3", lwd=2)
text(x = c(pts1-.3, pts2-.3, pts3-.3, pts4-.3, pts5-.3), y = c(sin(pts1), sin(pts2), sin(pts3), sin(pts4), sin(pts5)), c("A", "B", "C*", "C1", "C2"))

pts6 <- pts5 + sin(pts5)
pts7 <- pts4 - sin(pts4)*((pts4-pts2)/(sin(pts4)-sin(pts2)))
pts8 <- pts5 - sin(pts5)*((pts6-pts5)/(sin(pts6)-sin(pts5)))

x <- seq(0, 2*pi, .01)
plot(x, sin(x), "l", xlab = "x", ylab = "G(x)", main = "Iteration 2")
points(x = c(pts2, pts4, pts5, pts6, pts7, pts8), y = c(sin(pts2), sin(pts4), sin(pts5), sin(pts6), sin(pts7), sin(pts8)),  pch=19, cex=1.5)
abline(v = pi, lty=2)

segments(pts2, sin(pts2), pts4, sin(pts4), lty=3, col = "red", lwd=2)
segments(pts5, sin(pts5), pts6, sin(pts6), lty=2, col = "green3", lwd=2)
text(x = c(pts2-.3, pts4-.3, pts5-.3, pts6-.3, pts7+.3, pts8+.3), y = c(sin(pts2), sin(pts4), sin(pts5), sin(pts6), sin(pts7), sin(pts8)), c("B", "C1", "C2", "C**", "C3", "C4"))
dev.off()

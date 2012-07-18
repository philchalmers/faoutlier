library(faoutlier)

data(holzinger)
data(holzinger.outlier)
nfact <- 3

(MD <- robustMD(holzinger))
summary(MD)

(MD2 <- robustMD(holzinger.outlier))
summary(MD2)

(r1a <- LD(holzinger, nfact))
(r2a <- gCD(holzinger, nfact))
(r3a <- obs.resid(holzinger, nfact))
(r4a <- forward.search(holzinger, nfact))
plot(r1a)
plot(r2a)
plot(r3a)
plot(r4a)

(r1b <- LD(holzinger.outlier, nfact))
(r2b <- gCD(holzinger.outlier, nfact))
(r3b <- forward.search(holzinger.outlier, nfact))
(r4b <- obs.resid(holzinger.outlier, nfact))
plot(r1b)
plot(r2b)
plot(r3b)
plot(r4b)

###########

#one factor CFA with OpenMx
manifests <- colnames(holzinger)
latents <- c("G")
model <- mxModel("One Factor",
      type="RAM",
      manifestVars = manifests,
      latentVars = latents,
      mxPath(from=latents, to=manifests),
      mxPath(from=manifests, arrows=2),
      mxPath(from=latents, arrows=2,
            free=FALSE, values=1.0),
      mxData(cov(holzinger), type="cov", numObs=nrow(holzinger))
	  )			
summary(mxRun(model))

(r1 <- LD(holzinger, model))
(r2 <- gCD(holzinger, model))
(r3 <- obs.resid(holzinger, model))
(r4 <- forward.search(holzinger, model))
plot(r1)
plot(r2)
plot(r3)
plot(r4)


manifests <- colnames(holzinger)
latents <- c("G")
model <- mxModel("One Factor",
      type="RAM",
      manifestVars = manifests,
      latentVars = latents,
      mxPath(from=latents, to=manifests),
      mxPath(from=manifests, arrows=2),
      mxPath(from=latents, arrows=2,
            free=FALSE, values=1.0),
      mxData(cov(holzinger.outlier), type="cov", numObs=nrow(holzinger.outlier))
	  )			
summary(mxRun(model))


(r1 <- LD(holzinger.outlier, model))
(r2 <- gCD(holzinger.outlier, model))
(r3 <- obs.resid(holzinger.outlier, model))
(r4 <- forward.search(holzinger.outlier, model))
plot(r1)
plot(r2)
plot(r3)
plot(r4)

#three factor
manifests <- colnames(holzinger)
latents <- c("F1","F2","F3")
model3 <- mxModel("Three Factor",
      type="RAM",
      manifestVars = manifests,
      latentVars = latents,
      mxPath(from="F1", to=manifests[1:3]),
	  mxPath(from="F2", to=manifests[4:6]),
	  mxPath(from="F3", to=manifests[7:9]),
      mxPath(from=manifests, arrows=2),
      mxPath(from=latents, arrows=2,
            free=FALSE, values=1.0),
      mxData(cov(holzinger), type="cov", numObs=nrow(holzinger))
	  )			
summary(mxRun(model3))

(r1 <- LD(holzinger, model3))
(r2 <- gCD(holzinger, model3))
(r3 <- obs.resid(holzinger, model3))
(r4 <- forward.search(holzinger, model3))
plot(r1)
plot(r2)
plot(r3)
plot(r4)


#three factor with outlier
manifests <- colnames(holzinger)
latents <- c("F1","F2","F3")
model3 <- mxModel("Three Factor",
      type="RAM",
      manifestVars = manifests,
      latentVars = latents,
      mxPath(from="F1", to=manifests[1:3]),
	  mxPath(from="F2", to=manifests[4:6]),
	  mxPath(from="F3", to=manifests[7:9]),
      mxPath(from=manifests, arrows=2),
      mxPath(from=latents, arrows=2,
            free=FALSE, values=1.0),
      mxData(cov(holzinger.outlier), type="cov", numObs=nrow(holzinger.outlier))
	  )			
summary(mxRun(model3))

(r1 <- LD(holzinger.outlier, model3))
(r2 <- gCD(holzinger.outlier, model3))
(r3 <- obs.resid(holzinger.outlier, model3))
(r4 <- forward.search(holzinger.outlier, model3))
plot(r1)
plot(r2)
plot(r3)
plot(r4)

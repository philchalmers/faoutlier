context('obs.resid')

test_that('obs.resid run', {
    
    #Exploratory
    nfact <- 3
    ORresult <- obs.resid(holzinger, nfact)
    ORresult.outlier <- obs.resid(holzinger.outlier, nfact)
    expect_equal(ORresult$std_res[1:3], c(-0.01285074, -0.75147106, -1.42004924),
                 tolerance = 1e-5)
    expect_is(ORresult, 'obs.resid')
    expect_is(ORresult.outlier, 'obs.resid')
    expect_is(plot(ORresult), 'trellis')
    expect_is(plot(ORresult.outlier), 'trellis')    
    
    #-------------------------------------------------------------------    
    suppressMessages(model <- specifyModel(file='sem-model/sem-model.txt', quiet=TRUE))    
    ORresult <- obs.resid(holzinger, model)
    ORresult.outlier <- obs.resid(holzinger.outlier, model)
    expect_equal(ORresult$std_res[1:3], c(0.2548177, -0.5300287, -1.8518586),
                 tolerance = 1e-5)
    expect_is(ORresult, 'obs.resid')
    expect_is(ORresult.outlier, 'obs.resid')
    expect_is(plot(ORresult), 'trellis')
    expect_is(plot(ORresult.outlier), 'trellis')
    
    #-------------------------------------------------------------------    
    model <- 'F1 =~  Remndrs + SntComp + WrdMean
    F2 =~ MissNum + MxdArit + OddWrds
    F3 =~ Boots + Gloves + Hatchts'
    
    obs.resid2 <- obs.resid(holzinger, model, orthogonal=TRUE)
    obs.resid2.outlier <- obs.resid(holzinger.outlier, model, orthogonal=TRUE)
    expect_equal(ORresult$std_res[1:3], c(0.2548177, -0.5300287, -1.8518586),
                 tolerance = 1e-5)
    expect_is(obs.resid2, "obs.resid")
    expect_is(obs.resid2.outlier, "obs.resid")
    expect_is(plot(obs.resid2), 'trellis')
    expect_is(plot(obs.resid2.outlier), 'trellis')
})

test_that('gCD categorical', {
    set.seed(1234)
    cut <- rnorm(ncol(holzinger.outlier), 0 , .25)
    dat <- holzinger.outlier
    for(i in 1:length(cut))
        dat[,i] <- ifelse(holzinger.outlier[,i] < cut[i], 0, 1)
    
    dat <- as.data.frame(dat)
    model <- 'F1 =~  Remndrs + SntComp + WrdMean
    F2 =~ MissNum + MxdArit + OddWrds
    F3 =~ Boots + Gloves + Hatchts'
    
    obs.resid2 <- obs.resid(dat, model, orthogonal=TRUE, ordered=colnames(dat))
    expect_is(obs.resid2, "obs.resid")
    out <- print(obs.resid2)
    expect_equal(head(out), c(7.936152, 2.542361, 3.686347, 41.886093, 43.566606, 43.181056),
                 tolerance = 1e-4)
})
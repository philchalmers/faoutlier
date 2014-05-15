context('gCD')

test_that('gCD run', {
    
    #Exploratory
    nfact <- 3
    gCDresult <- gCD(holzinger, nfact)
    gCDresult.outlier <- gCD(holzinger.outlier, nfact)
    expect_is(gCDresult, 'gCD')
    expect_is(plot(gCDresult), 'trellis')
    expect_is(gCDresult.outlier, 'gCD')
    expect_is(plot(gCDresult.outlier), 'trellis')
    
    #-------------------------------------------------------------------    
    suppressMessages(model <- specifyModel(file='sem-model/sem-model.txt', quiet=TRUE))
    gCDresult <- gCD(holzinger, model)
    gCDresult.outlier <- suppressWarnings(gCD(holzinger.outlier, model))
    expect_is(gCDresult, 'gCD')
    expect_is(plot(gCDresult), 'trellis')
    expect_is(gCDresult.outlier, 'gCD')
    expect_is(plot(gCDresult.outlier), 'trellis')
    
    #-------------------------------------------------------------------
    #Confirmatory with lavaan
    model <- 'F1 =~  Remndrs + SntComp + WrdMean
    F2 =~ MissNum + MxdArit + OddWrds
    F3 =~ Boots + Gloves + Hatchts'
    
    gCDresult <- gCD(holzinger, model, orthogonal=TRUE)
    gCDresult.outlier <- gCD(holzinger.outlier, model, orthogonal=TRUE)
    expect_equal(head(gCDresult.outlier$gCD), c(0.0104, 0.0000, 0.0001, 0.0001, 0.0000, 0.0002),
                 tolerance = 1e-2)
    expect_is(gCDresult, 'gCD')
    expect_is(plot(gCDresult), 'trellis')
    expect_is(gCDresult.outlier, 'gCD')
    expect_is(plot(gCDresult.outlier), 'trellis')
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
    
    gCDresult <- gCD(dat, model, orthogonal=TRUE, ordered=colnames(dat))
    expect_is(gCDresult, 'gCD')
    expect_equal(head(gCDresult$gCD), c(0.0007, 0.0005, 0.0069, 0.0005, 0.0012, 0.0038),
                 tolerance = 1e-2)
})


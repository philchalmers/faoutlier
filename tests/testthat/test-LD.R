context('LD')

test_that('LD run', {    
    
    #Exploratory
    nfact <- 3
    LDresult <- LD(holzinger, nfact)    
    LDresult.outlier <- LD(holzinger.outlier, nfact)
    expect_equal(as.numeric(LDresult[1:3]), c(-0.7431551, -0.3786938, -1.2646524),
                 tolerance=1e-5)
    expect_is(LDresult, 'LD')
    expect_is(LDresult.outlier, 'LD')
    expect_is(plot(LDresult), 'trellis')
    expect_is(plot(LDresult.outlier), 'trellis')
    
    #-------------------------------------------------------------------    
    suppressMessages(model <- specifyModel(file='sem-model/sem-model.txt', quiet=TRUE))
    LDresult <- LD(holzinger, model)
    LDresult.outlier <- LD(holzinger.outlier, model)
    expect_equal(as.numeric(LDresult[1:3]), c(-4.945440, -1.943843, -3.330193),
                 tolerance=1e-5)
    expect_is(LDresult, 'LD')
    expect_is(LDresult.outlier, 'LD')
    expect_is(plot(LDresult), 'trellis')
    expect_is(plot(LDresult.outlier), 'trellis')
    
    #-------------------------------------------------------------------
    #Confirmatory with lavaan
    model <- 'F1 =~  Remndrs + SntComp + WrdMean
    F2 =~ MissNum + MxdArit + OddWrds
    F3 =~ Boots + Gloves + Hatchts'
    
    LDresult <- LD(holzinger, model, orthogonal=TRUE)
    expect_equal(as.numeric(LDresult[1:3]), c(-4.984276, -1.952360, -3.352713),
                 tolerance=1e-5)
    LDresult <- LD(holzinger.outlier, model, orthogonal=TRUE)
    expect_is(LDresult, 'LD')
    expect_is(LDresult.outlier, 'LD')
    expect_is(plot(LDresult), 'trellis')
    expect_is(plot(LDresult.outlier), 'trellis')
})

test_that('LD categorical', {
    set.seed(1234)
    cut <- rnorm(ncol(holzinger.outlier), 0 , .25)
    dat <- holzinger.outlier
    for(i in 1:length(cut))
        dat[,i] <- ifelse(holzinger.outlier[,i] < cut[i], 0, 1)
    
    dat <- as.data.frame(dat)
    model <- 'F1 =~  Remndrs + SntComp + WrdMean
    F2 =~ MissNum + MxdArit + OddWrds
    F3 =~ Boots + Gloves + Hatchts'
    
    LDresult <- LD(dat, model, orthogonal=TRUE, ordered=colnames(dat))
    expect_is(LDresult, 'LD')
    expect_equal(as.numeric(head(LDresult)), 
                 c(2.392693, -8.217107, -2.637149, -17.303106, 4.298777, -11.997538),
                 tolerance = 1e-3)
})
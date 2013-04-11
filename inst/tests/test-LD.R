context('LD')

test_that('LD run', {    
    
    #Exploratory
    nfact <- 3
    LDresult <- LD(holzinger, nfact)    
    (LDresult.outlier <- LD(holzinger.outlier, nfact))
    expect_is(LDresult, 'LD')
    expect_is(LDresult.outlier, 'LD')
    expect_is(plot(LDresult), 'trellis')
    expect_is(plot(LDresult.outlier), 'trellis')
    
    #-------------------------------------------------------------------    
    suppressMessages(model <- specifyModel(file='sem-model.txt', quiet=TRUE))
    LDresult <- LD(holzinger, model)
    LDresult.outlier <- LD(holzinger.outlier, model)
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
    LDresult <- LD(holzinger.outlier, model, orthogonal=TRUE)
    expect_is(LDresult, 'LD')
    expect_is(LDresult.outlier, 'LD')
    expect_is(plot(LDresult), 'trellis')
    expect_is(plot(LDresult.outlier), 'trellis')
})
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
    suppressMessages(model <- specifyModel(file='sem-model.txt', quiet=TRUE))
    gCDresult <- gCD(holzinger, model)
    gCDresult.outlier <- gCD(holzinger.outlier, model)
    expect_is(gCDresult, 'gCD')
    expect_is(plot(gCDresult), 'trellis')
    expect_is(gCDresult.outlier, 'gCD')
    expect_is(plot(gCDresult.outlier), 'trellis')
    
    #-------------------------------------------------------------------
    #Confirmatory with lavaan
    model <- 'F1 =~  Remndrs + SntComp + WrdMean
    F2 =~ MissNum + MissNum + OddWrds
    F3 =~ Boots + Gloves + Hatchts'
    
    gCDresult <- gCD(holzinger, model)
    gCDresult.outlier <- gCD(holzinger.outlier, model)
    expect_is(gCDresult, 'gCD')
    expect_is(plot(gCDresult), 'trellis')
    expect_is(gCDresult.outlier, 'gCD')
    expect_is(plot(gCDresult.outlier), 'trellis')
})
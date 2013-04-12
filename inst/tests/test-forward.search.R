context('forward.search')

test_that('forward.search run', {    
    
    #Exploratory
    nfact <- 3
    FS <- forward.search(holzinger, nfact, print.messages = FALSE)    
    expect_is(FS, 'forward.search')
    expect_is(plot(FS), 'trellis')
    
    #-------------------------------------------------------------------    
    suppressMessages(model <- specifyModel(file='sem-model.txt', quiet=TRUE))    
    FS.outlier <- forward.search(holzinger.outlier, model, print.messages = FALSE)
    expect_is(FS.outlier, 'forward.search')
    expect_is(plot(FS.outlier), 'trellis')
    
    #---- lavaan
    model <- 'F1 =~  Remndrs + SntComp + WrdMean
    F2 =~ MissNum + MxdArit + OddWrds
    F3 =~ Boots + Gloves + Hatchts'
    FS.outlier <- forward.search(holzinger.outlier, model, print.messages = FALSE)
    expect_is(FS.outlier, 'forward.search')
    expect_is(plot(FS.outlier), 'trellis')
    
})
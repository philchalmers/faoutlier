context('gCD')

test_that('gCD run', {

    #Exploratory
    nfact <- 3
    gCDresult <- gCD(holzinger, nfact, progress = FALSE, vcov_drop = TRUE)
    gCDresult.outlier <- gCD(holzinger.outlier, nfact, progress = FALSE)
    expect_is(gCDresult, 'gCD')
    expect_is(plot(gCDresult), 'trellis')
    expect_is(gCDresult.outlier, 'gCD')
    expect_is(plot(gCDresult.outlier), 'trellis')
    expect_equal(as.numeric(gCDresult$gCD[1:3]), c(0.0020993686, 0.0008613305, 0.0013162123), tolerance=1e-5)
    gCDresult <- gCD(holzinger, nfact, progress = FALSE)
    expect_equal(as.numeric(gCDresult$gCD[1:3]), c(0.0022232819, 0.0008943396, 0.0012875196), tolerance=1e-5)

    #-------------------------------------------------------------------
    suppressMessages(model <- sem::specifyModel(file='sem-model/sem-model.txt', quiet=TRUE))
    gCDresult <- suppressWarnings(gCD(holzinger, model, progress = FALSE, vcov_drop = TRUE))
    gCDresult.outlier <- suppressWarnings(gCD(holzinger.outlier, model, progress = FALSE, , vcov_drop = TRUE))
    expect_is(gCDresult, 'gCD')
    expect_is(plot(gCDresult), 'trellis')
    expect_is(gCDresult.outlier, 'gCD')
    expect_is(plot(gCDresult.outlier), 'trellis')
    expect_equal(as.numeric(gCDresult$gCD[1:3]), c(0.14419229, 0.04918079, 0.42079870), tolerance=1e-5)

    #-------------------------------------------------------------------
    #Confirmatory with lavaan
    model <- 'F1 =~  Remndrs + SntComp + WrdMean
    F2 =~ MissNum + MxdArit + OddWrds
    F3 =~ Boots + Gloves + Hatchts'

    gCDresult <- gCD(holzinger, model, orthogonal=TRUE, progress = FALSE, vcov_drop = TRUE)
    gCDresult.outlier <- gCD(holzinger.outlier, model, orthogonal=TRUE, progress = FALSE, vcov_drop = TRUE)
    expect_equal(as.numeric(gCDresult.outlier$gCD[1:3]), c(36.67371343, 0.05034515, 0.34999495),
                 tolerance = 1e-5)
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

    gCDresult <- suppressWarnings(gCD(dat, model, orthogonal=TRUE, progress = FALSE,
                                      vcov_drop = TRUE, ordered=colnames(dat)))
    expect_is(gCDresult, 'gCD')
    expect_equal(as.numeric(head(gCDresult$gCD)), c(0.45389878, 0.11746119, 0.15885323, 0.07169155, 0.30643962, 0.15728745),
                 tolerance = 1e-2)

    #model <- mirt::mirt.model('F1 = 1-3
    #                         F2 = 4-6
    #                         F3 = 7-9')
    #gCDresult2 <- suppressWarnings(gCD(dat, model))

})


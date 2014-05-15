context('robustMD')

test_that('robustMD run', {    
    output <- robustMD(holzinger)
    expect_is(output, 'robmah')
    expect_is(plot(output), 'trellis')
    expect_is(plot(output, type = 'qqplot'), 'trellis')
})
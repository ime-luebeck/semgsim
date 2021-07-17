context("sample_MU_contributions")

test_that("sample_randomly_shifted_time_response works - simple test case", {

    res <- sample_randomly_shifted_time_response(
        firing_idx = 10,
        TF_fft_vals = c(2, 5+1i, 3-1i, 2, 0, 2, 3+1i, 5-1i),
        sampling = create.sampling(Fs = 1000, NFFT = 8),
        n_samples = 40)
    
    ## The resulting sequence is 0 everywhere except on the support of the
    ## firing response.
    expect_true(all(res[-(10:17)] == 0))
    expect_equal(length(res), 40)

    res2 <- sample_randomly_shifted_time_response(
        firing_idx = 10,
        TF_fft_vals = c(2, 5+1i, 3-1i, 2, 1i, -1i, 2, 3+1i, 5-1i),
        sampling = create.sampling(Fs = 1000, NFFT = 9),
        n_samples = 40)
    
    ## The resulting sequence is 0 everywhere except on the support of the
    ## firing response.
    expect_true(all(res2[-(10:18)] == 0))
    expect_equal(length(res2), 40)
})

test_that("sample_randomly_shifted_time_response throws errors correctly", {

    ## We expect an error here since the TF_fft_vals are not hermitian, which is
    ## a prerequisite.
    expect_error(
        sample_randomly_shifted_time_response(
            firing_idx = 10,
            TF_fft_vals = c(1, 2, 3+i, 4, 0, 4, 3, 2),
            sampling =
                create.sampling(
                    Fs = 1000,
                    NFFT = 8),
            n_samples = 40))

    ## Here, the fft is hermitian (since resulting from a real signal) but has a
    ## non-zero coefficient at Fs/2, which leads to problems with the time shift
    ## operation. Thus we expect the function to throw an error.
    expect_error(
        sample_randomly_shifted_time_response(
            firing_idx = 10,
            TF_fft_vals = fft(1:10),
            sampling =
                create.sampling(Fs = 1000,
                                NFFT = 10),
            n_samples = 40))
    
})

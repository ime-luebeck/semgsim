context("sample_MU_contribs_to_potential")

NFFT <- 100
Fs <- 100
sampling <- create_sampling(Fs = Fs, NFFT = NFFT)

firing_responses_fftvals <- fft(rnorm(100))
firing_instants <- sample(c(T,F), size = Fs * 10, replace = TRUE)

test("sequential and parallel implementations yield the same result", {

    ## Can't due to randomness being involved. Would require further
    ## infrastructure.
})

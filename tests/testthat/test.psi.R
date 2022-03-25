test_that("psi seems implemented correctly", {
  #' This should be d/dz V(-z)
  library(numDeriv)
  z_m <- seq(-0.03, 0.03, 1/1000)
  Vminus <- function(z) IAP_Rosenfalck(-z)
  suppressWarnings(Vminus_diff <- grad(Vminus, z_m))
  psi <- psi_Rosenfalck(z_m)
  expect_equal(psi, Vminus_diff, 1e-5)
})

test_that("psi_trafo seems implemented correctly", {
  #' This should be d/dz V(-z)
  
  Fz <- 4096
  NFFT <- 256
  sampling <- semgsim:::create_sampling(Fz, NFFT)
  z <- seq(-0.0625 + 1/Fz, 0, 1/Fz)
  psi_vals <- psi_Rosenfalck(z)
  psi_vals_fft <- fft(psi_vals) / Fz
  psi_vals_recovered <- Re(fft(psi_vals_fft, inverse = TRUE) * Fz / NFFT)
  # this does not test psi_trafo but just whether the forward-backward FFT scheme is right
  expect_equal(psi_vals, psi_vals_recovered)
  
  psi_trafo_vals_fft <- psi_Rosenfalck_transformed(sampling$fft_freqs)
  # This test does not pass, although the abs(fft) curves are very similar.
  # However, the Re and Im differ significantly. Is this a problem?
  # If yes, which of the two is right?
  expect_equal(psi_vals_fft, psi_trafo_vals_fft)
})

test_that("psi_trafo seems implemented correctly #1", {
  #' This should be d/dz V(-z)
  
  Fs <- 1024
  NFFT <- 2048
  sampling <- semgsim:::create_sampling(Fs, NFFT)
  
  calc_TF_point <- function(ang_freq) {
    real_part_1 <- try(integrate(function(z) Re(psi_Rosenfalck(z) * exp(-1i * ang_freq * z)), -Inf, 0)$value, silent = TRUE)
    real_part_2 <- try(integrate(function(z) Re(psi_Rosenfalck(z) * exp(-1i * ang_freq * z)), 0, Inf)$value, silent = TRUE)
    if (class(real_part_1) == 'try-error' || class(real_part_2) == 'try-error') {
      real_part <- 0
    } else
      real_part <- real_part_1 + real_part_2
    
    imaginary_part_1 <- try(integrate(function(z) Im(psi_Rosenfalck(z) * exp(-1i * ang_freq * z)), -Inf, 0)$value, silent = TRUE)
    imaginary_part_2 <- try(integrate(function(z) Im(psi_Rosenfalck(z) * exp(-1i * ang_freq * z)), 0, Inf)$value, silent = TRUE)
    if (class(imaginary_part_1) == 'try-error' || class(imaginary_part_2) == 'try-error') {
      imaginary_part <- 0
    } else
      imaginary_part <- imaginary_part_1 + imaginary_part_2
    
    return(real_part + 1i * imaginary_part)
  }
  calc_TF <- Vectorize(calc_TF_point, 'ang_freq')
  
  psi_trafo_vals_analytical <- calc_TF(sampling$freqs_to_calc * 2 * pi)
  psi_trafo_vals <- psi_Rosenfalck_transformed(sampling$freqs_to_calc)

  expect_equal(psi_trafo_vals_analytical, psi_trafo_vals)
})



test_that("psi_trafo seems implemented correctly #2", {
  #' This should be d/dz V(-z)
  
  Fz <- 1024
  NFFT <- 2048
  sampling <- semgsim:::create_sampling(Fz, NFFT)
  z = seq(-1, 1, 1/Fz)
  # We only need to calculate half of the FFT values, because for a real signal, the Fourier spectrum is Hermitian
  psi_trafo_vals <- psi_Rosenfalck_transformed(sampling$freqs_to_calc)
  
  # Obtain the remaining samples by mirroring
  fftvals <- semgsim:::mirror_hermitian_fft(psi_trafo_vals, sampling)
  sig_reconstruction <- fft(fftvals, inverse=TRUE) * Fz / NFFT
  orig_signal = psi_Rosenfalck(z)
  expect_equal(sig_reconstruction, psi_trafo_vals)
})
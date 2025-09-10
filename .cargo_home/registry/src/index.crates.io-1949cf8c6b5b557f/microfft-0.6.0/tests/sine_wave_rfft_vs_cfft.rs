use microfft::Complex32;
use std::convert::TryInto;
use std::f32::consts::PI;

/// This test produces 1024 samples of a 129 Hertz sine wave with a
/// sampling rate of 44.1kHz and applies real and complex FFT on it.
#[test]
fn test_real_vs_complex_fft_with_real_sine_wave() {
    const SAMPLE_COUNT: usize = 1024;
    let sampling_rate = 44100.0;
    let frequency_resolution = sampling_rate / SAMPLE_COUNT as f32;
    // Frequency of test signal in Hertz
    let test_frequency = 129.0;

    let sine_wave_fn = sine_wave(test_frequency);

    let mut samples_real: [f32; SAMPLE_COUNT] =
        get_sine_wave_samples::<SAMPLE_COUNT>(&sine_wave_fn, sampling_rate);
    let mut samples_complex: [Complex32; SAMPLE_COUNT] = samples_real
        .iter()
        .map(|f| Complex32::new(*f, 0.0))
        .collect::<Vec<Complex32>>()
        .try_into()
        .unwrap();

    // frequency_resolution == sampling_rate / SAMPLE_COUNT == 44100 / 1024 == ~43 Hertz
    //
    // frequency == index_fft_result * frequency_resolution
    // 129 Hertz == index_fft_result * 43 Hertz => index_fft_result == 3

    let rfft_res = microfft::real::rfft_1024(&mut samples_real);
    let cfft_res = microfft::complex::cfft_1024(&mut samples_complex);

    let rfft_res_f32: [f32; SAMPLE_COUNT / 2] = rfft_res
        .iter()
        .map(|c| c.norm() / SAMPLE_COUNT as f32)
        .collect::<Vec<f32>>()
        .try_into()
        .unwrap();
    let cfft_res_f32: [f32; SAMPLE_COUNT] = cfft_res
        .iter()
        .map(|c| c.norm() / SAMPLE_COUNT as f32)
        .collect::<Vec<f32>>()
        .try_into()
        .unwrap();

    // print spectrum 0 to 172 Hertz
    // WE EXPECT A MAXIMUM/PEAK AT 129 Hertz!
    for i in 0..5 {
        let fr = i as f32 * frequency_resolution;
        println!("RealFFT   : {}Hz => {}", fr, rfft_res_f32[i]);
        println!("ComplexFFT: {}Hz => {}", fr, cfft_res_f32[i]);
    }

    // Assert clearly detectable maximum/peak
    assert!(
        rfft_res_f32[2] < 2.0 && rfft_res_f32[3] > 500.0 && rfft_res_f32[4] < 3.0,
        "The peak must be exactly at 129 Hertz!"
    );

    // assert real FFT and complex FFT are equal
    for i in 0..5 {
        assert!(
            (rfft_res_f32[i] - cfft_res_f32[i]).abs() <= 0.01,
            "Outputs of real and complex FFTs must be equal"
        );
    }
}

/// Helper function for test [`test_real_vs_complex_fft_with_real_sine_wave`].
/// Generates a sine wave function of a given frequency.
fn sine_wave(frequency: f32) -> impl Fn(f32) -> f32 {
    move |t| (t * frequency * 2.0 * PI).sin()
}

/// Helper function for test [`test_real_vs_complex_fft_with_real_sine_wave`].
/// Generates N samples using a sine wave function. Or in other words:
/// This produces audio data with a sine wave of a given frequency over
/// the time. The time is implicitly determined by the number of samples
/// and the sampling rate.
fn get_sine_wave_samples<const N: usize>(
    sine_wave_fn: &dyn Fn(f32) -> f32,
    sampling_rate: f32,
) -> [f32; N] {
    let mut samples = [0.0; N];
    for i in 0..N {
        let t = (1.0 / sampling_rate as f32) * i as f32;
        samples[i] = sine_wave_fn(t as f32) * 1000.0;
    }
    samples
}

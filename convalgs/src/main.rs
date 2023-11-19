pub mod convalgs;

use convalgs::{Convolve, ConvolutionMethod};
use rand::{thread_rng, Rng};
use std::time::Instant;


fn main() {

    let mut x: Vec<f32> = vec![0.0; 44100];
    let mut h: Vec<f32> = vec![0.0; 3000];

    let mut rng = thread_rng();
    let x = x.iter_mut().map(|_| rng.gen_range(-1.0..1.0)).collect();
    let h = h.iter_mut().map(|_| rng.gen_range(-1.0..1.0)).collect();

    let conv = Convolve::new(&x, &h);
    
    get_process_time(&conv, ConvolutionMethod::InputSide, "INPSIDE");
    get_process_time(&conv, ConvolutionMethod::OutputSide, "OUTSIDE");
    get_process_time(&conv, ConvolutionMethod::Karatsuba, "KARAT");
    get_process_time(&conv, ConvolutionMethod::Fft, "FFT");
    get_process_time(&conv, ConvolutionMethod::OlaFft(4096), "OLA FFT");

}

fn get_process_time(conv_planner: &Convolve<f32>, conv_method: ConvolutionMethod, label: &str) {
    let now = Instant::now();
    let _ = conv_planner.convolve(conv_method);
    let elapsed_time = now.elapsed().as_secs_f64();
    println!("{} elapsed time: {}", label, elapsed_time);
}
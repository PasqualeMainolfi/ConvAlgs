pub mod convalgs;
use convalgs::{Convolve, ConvolutionMethod};


fn main() {

    let x: Vec<f32> = [0.1, -3.2, 1.5, 7.0, 5.7].to_vec();
    let h: Vec<f32> = [0.7, -3.4, 1.0].to_vec();

    let conv = Convolve::new(&x, &h);
    let inpside = conv.convolve(ConvolutionMethod::InputSide);
    let outside = conv.convolve(ConvolutionMethod::OutputSide);
    let karat = conv.convolve(ConvolutionMethod::Karatsuba);
    let fft = conv.convolve(ConvolutionMethod::Fft);

    println!("INPSIDE: {:?}", inpside);
    println!("OUTSIDE: {:?}", outside);
    println!("KARAT: {:?}", karat);
    println!("FFT: {:?}", fft);

}

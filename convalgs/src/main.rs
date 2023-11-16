pub mod convalgs;
use convalgs::{Convolve, ConvolutionMethod};


fn main() {
    let x: Vec<f32> = [0.1, -3.2, 1.5, 7.0, 5.7].to_vec();
    let h: Vec<f32> = [0.1, -3.2, 1.5].to_vec();

    let conv = Convolve::new(&x, &h);
    let convolution = conv.convolve(ConvolutionMethod::OutputSide);
    println!("{:?}", convolution);
}


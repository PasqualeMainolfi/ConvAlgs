use num_traits::Float;


pub enum ConvolutionMethod {
    OutputSide,
    InputSide,
    Karatsuba,
    Fft,
    OlaFft,
    Ntt
}

pub struct Convolve<'a, T: Float>
where
    T: Float 
{
    x: &'a Vec<T>,
    h: &'a Vec<T>,
}

impl<'a, T> Convolve<'a, T>
where 
    T: Float 
{
    pub fn new(x: &'a Vec<T>, h: &'a Vec<T>) -> Self {
        Self { x, h, }
    }

    pub fn convolve(&self, mode: ConvolutionMethod) -> Vec<T> {
        let mut y: Vec<T> = Vec::new();
        match mode {
            ConvolutionMethod::InputSide => { todo!() },
            ConvolutionMethod::OutputSide => { todo!() },
            ConvolutionMethod::Karatsuba => { todo!() },
            ConvolutionMethod::Fft => { todo!() },
            ConvolutionMethod::OlaFft => { todo!() },
            ConvolutionMethod::Ntt => { todo!() },
        }
    }



}

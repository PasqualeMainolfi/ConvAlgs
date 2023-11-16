use num_traits::Float;
use std::ops::{Add, Mul, Div, Sub, AddAssign};


pub enum ConvolutionMethod {
    OutputSide,
    InputSide,
    Karatsuba,
    Fft,
    OlaFft,
    Ntt
}

pub struct Convolve<'a, T: Float> {
    x: &'a Vec<T>,
    h: &'a Vec<T>,
}

impl<'a, T> Convolve<'a, T> 
where
    T: Float
{
    pub fn new(x: &'a Vec<T>, h: &'a Vec<T>) -> Self {
        Self { x, h }
    }

    pub fn convolve(&self, mode: ConvolutionMethod) -> Vec<T> 
    where
        T: Float + Add + Mul + Div + Sub + AddAssign
    {
        let xlen = self.x.len();
        let hlen = self.h.len();
        let ylen = xlen + hlen - 1;
        let mut y: Vec<T> = vec![T::zero(); ylen];

        match mode {
            ConvolutionMethod::InputSide => {
                for k in 0..hlen {
                    for i in 0..xlen {
                        y[k + i] += self.h[k] * self.x[i]
                    }
                }
            },
            ConvolutionMethod::OutputSide => {
                for (n, value) in y.iter_mut().enumerate() {
                    let lower = 0.max(n as isize - xlen as isize + 1) as usize;
                    let upper = n.min(hlen - 1);
                    for k in lower..upper {
                        *value += self.h[k] * self.x[n - k];
                    }
                }
            },
            ConvolutionMethod::Karatsuba => { todo!() },
            ConvolutionMethod::Fft => { todo!() },
            ConvolutionMethod::OlaFft => { todo!() },
            ConvolutionMethod::Ntt => { todo!() },
        }
        y
    }

}

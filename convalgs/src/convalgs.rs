#![allow(clippy::only_used_in_recursion)]

use num_traits::Float;
use std::ops::{Add, Mul, AddAssign};
use std::fmt::Debug;
use realfft::{RealFftPlanner, FftNum};
use rustfft::num_complex::Complex;


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
    T: Debug + Float + AddAssign + Add + Mul + FftNum +
{
    pub fn new(x: &'a Vec<T>, h: &'a Vec<T>) -> Self {
        Self { x, h }
    }

    pub fn convolve(&self, mode: ConvolutionMethod) -> Vec<T> {
        let xlen = self.x.len();
        let hlen = self.h.len();
        let ylen = xlen + hlen - 1;
        let mut y: Vec<T> = vec![T::zero(); ylen];

        let pow_size = 1 << (ylen as f32).log2().ceil() as usize;
        let mut xpad = vec![T::zero(); pow_size];
        let mut hpad = vec![T::zero(); pow_size];
        xpad[..xlen].copy_from_slice(&self.x[..]);
        hpad[..hlen].copy_from_slice(&self.h[..]);

        match mode {
            ConvolutionMethod::InputSide => self._input_side_helper(self.x, self.h, &mut y),
            ConvolutionMethod::OutputSide => self._output_side_helper(self.x, self.h, &mut y),
            ConvolutionMethod::Karatsuba => {
                let _y = self._karatsuba_helper(&xpad, &hpad);
                y[..].copy_from_slice(&_y[..ylen]);
            },
            ConvolutionMethod::Fft => self._fft_helper(&mut xpad, &mut hpad, &mut y),
            ConvolutionMethod::OlaFft => { todo!() },
            ConvolutionMethod::Ntt => { todo!() },
        }
        y
    }

    fn _input_side_helper(&self, x: &[T], h: &[T], buffer: &mut [T]) {
        for k in 0..h.len() {
            for i in 0..x.len() {
                buffer[k + i] += h[k] * x[i]
            }
        }
    }
    
    fn _output_side_helper(&self, x: &[T], h: &[T], buffer: &mut [T]) {
        for (n, value) in buffer.iter_mut().enumerate() {
            let lower = 0.max(n as isize - x.len() as isize + 1) as usize;
            let upper = n.min(h.len() - 1);
            for k in lower..=upper {
                *value += h[k] * x[n - k];
            }
        }
    }

    fn _karatsuba_helper(&self, x: &Vec<T>, h: &Vec<T>) -> Vec<T> {
        let xlen = x.len();
        let hlen = h.len();

        let n = xlen.max(hlen);

        if n == 1 { 
            let y = [x[0] * h[0]].to_vec();
            return y
        }

        let m = n / 2;
        let (x0, x1) = x.split_at(m);
        let (h0, h1) = h.split_at(m);

        let z0 = self._karatsuba_helper(&x0.to_vec(), &h0.to_vec());
        let z2 = self._karatsuba_helper(&x1.to_vec(), &h1.to_vec());

        let xsum = x0.iter().zip(x1.iter()).map(|(&a, &b)| a + b).collect::<Vec<_>>();
        let hsum = h0.iter().zip(h1.iter()).map(|(&a, &b)| a + b).collect::<Vec<_>>();

        let z1 = self._karatsuba_helper(&xsum, &hsum);
        let z1 = z1.iter().zip(z2.iter().zip(z0.iter())).map(|(&a, (&b, &c))| a - b - c).collect::<Vec<_>>();
        
        let mut y = vec![T::zero(); 2 * n - 1];
        for i in 0..z0.len() {
            y[i] += z0[i]
        }
        for i in 0..z1.len() {
            y[i + m] += z1[i]
        }
        for i in 0..z2.len() {
            y[i + 2 * m] += z2[i]
        }
        
        y

    }

    fn _fft_helper(&self, x: &mut [T], h: &mut [T], buffer: &mut [T]) {
        let ylen = x.len();
        
        let mut xplanner = RealFftPlanner::<T>::new();
        let mut hplanner = RealFftPlanner::<T>::new();
        
        let xfft = xplanner.plan_fft_forward(ylen);
        let hfft = hplanner.plan_fft_forward(ylen);

        let mut xspectrum = xfft.make_output_vec();
        let mut hspectrum = hfft.make_output_vec();

        xfft.process(x, &mut xspectrum).unwrap();
        hfft.process(h, &mut hspectrum).unwrap();

        let mut fft_prod = xspectrum.iter().zip(hspectrum.iter()).map(|(&a, &b)| a * b).collect::<Vec<Complex<T>>>();

        let ifft = xplanner.plan_fft_inverse(ylen);
        let mut ifft_time = ifft.make_output_vec();
        ifft.process(&mut fft_prod, &mut ifft_time).unwrap();

        for (i, value) in buffer.iter_mut().enumerate() {
            *value = ifft_time[i] / T::from(ylen).unwrap();
        }

    }

}

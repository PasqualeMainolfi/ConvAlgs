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
    OlaFft(usize),
}

pub struct Convolve<'a, T: Float> {
    x: &'a Vec<T>,
    h: &'a Vec<T>,
}

impl<'a, T> Convolve<'a, T> 
where
    T: Debug + Float + AddAssign + Add + Mul + FftNum
{
    pub fn new(x: &'a Vec<T>, h: &'a Vec<T>) -> Self {
        Self { x, h }
    }

    pub fn convolve(&self, mode: ConvolutionMethod) -> Vec<T> {
        let xlen = self.x.len();
        let hlen = self.h.len();
        let ylen = xlen + hlen - 1;
        let mut y: Vec<T> = vec![T::zero(); ylen];

        match mode {
            ConvolutionMethod::InputSide => self._input_side_helper(self.x, self.h, &mut y),
            ConvolutionMethod::OutputSide => self._output_side_helper(self.x, self.h, &mut y),
            ConvolutionMethod::Karatsuba => {
                let pow_size = 1 << (ylen as f32).log2().ceil() as usize;
                let mut xpad = vec![T::zero(); pow_size];
                let mut hpad = vec![T::zero(); pow_size];
                xpad[..xlen].copy_from_slice(&self.x[..]);
                hpad[..hlen].copy_from_slice(&self.h[..]);
                let _y = self._karatsuba_helper(&xpad, &hpad);
                y[..].copy_from_slice(&_y[..ylen]);
            },
            ConvolutionMethod::Fft => self._fft_helper(self.x, self.h, &mut y),
            ConvolutionMethod::OlaFft(value) => self._olafft_helper(&mut y, value),
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

    fn _fft_helper(&self, x: &[T], h: &[T], buffer: &mut [T]) {
        let xlen = x.len();
        let hlen = h.len();
        let ylen = xlen + hlen - 1;
        let conv_len = 1 << (ylen as f32).log2().ceil() as usize;

        let mut xpad = vec![T::zero(); conv_len];
        let mut hpad = vec![T::zero(); conv_len];

        xpad[..xlen].copy_from_slice(&x[..xlen]);
        hpad[..hlen].copy_from_slice(&h[..hlen]);

        let mut xplanner = RealFftPlanner::<T>::new();
        let mut hplanner = RealFftPlanner::<T>::new();
        
        let xfft = xplanner.plan_fft_forward(conv_len);
        let hfft = hplanner.plan_fft_forward(conv_len);

        let mut xspectrum = xfft.make_output_vec();
        let mut hspectrum = hfft.make_output_vec();

        xfft.process(&mut xpad, &mut xspectrum).unwrap();
        hfft.process(&mut hpad, &mut hspectrum).unwrap();

        let mut fft_prod = xspectrum.iter().zip(hspectrum.iter()).map(|(&a, &b)| a * b).collect::<Vec<Complex<T>>>();

        let ifft = xplanner.plan_fft_inverse(conv_len);
        let mut ifft_time = ifft.make_output_vec();
        ifft.process(&mut fft_prod, &mut ifft_time).unwrap();

        for (i, value) in buffer.iter_mut().enumerate() {
            *value = ifft_time[i] / T::from(conv_len).unwrap(); // return the buffer with len n + m - 1
        }

    }

    fn _olafft_helper(&self, buffer: &mut [T], frame_size: usize) {
        let xlen = self.x.len();
        let hlen = self.h.len();
        
        let nframes: usize = (xlen as f32 / frame_size as f32).ceil() as usize; // get number of total frames
        let xlen_new = nframes * frame_size; // adjust x len = nframes * frames size

        let len_fft_buffer = frame_size + hlen - 1; // len of fft buffer -> frame size + kernel size - 1
        let ylen = nframes * frame_size + hlen - 1; // len of result vector

        let mut xnewlen = vec![T::zero(); xlen_new];
        xnewlen[..xlen].copy_from_slice(&self.x[..xlen]);

        let mut fft_buffer: Vec<T> = vec![T::zero(); len_fft_buffer]; // fft buffer
        let mut y = vec![T::zero(); ylen]; // result vector
        
        for i in 0..nframes {
            let frame_start = i * frame_size;
            let frame_end = (i + 1) * frame_size;
            self._fft_helper(&xnewlen[frame_start..frame_end], self.h, &mut fft_buffer);

            for j in 0.. len_fft_buffer {
                y[frame_start + j] += fft_buffer[j]; // rebuild from ola with len nframes * frame size + hsize - 1
            }
        }

        for (i, value) in buffer.iter_mut().enumerate() {
            *value = y[i]; // return the buffer with len n + m - 1
        }

    }

}

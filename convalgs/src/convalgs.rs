#![allow(clippy::only_used_in_recursion)]

use num_traits::Float;
use std::ops::{Add, Mul, AddAssign};
use std::fmt::Debug;
use realfft::{RealFftPlanner, FftNum};
use rustfft::num_complex::Complex;

#[derive(Debug)]
pub enum ConvolutionMethod {
    OutputSide,
    InputSide,
    Fft,
    OlaFft(usize),
}

pub struct Convolve<'a, T> 
where
    T: Float
{
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
            ConvolutionMethod::Fft => {
                let pow_size = 1 << ((ylen as f32).log2() + 1.0) as usize;
                let mut xpad = vec![T::zero(); pow_size];
                let mut hpad = vec![T::zero(); pow_size];
                xpad[..xlen].copy_from_slice(&self.x[..]);
                hpad[..hlen].copy_from_slice(&self.h[..]);
                self._fft_helper(&mut xpad, &mut hpad, &mut y)
            },
            ConvolutionMethod::OlaFft(value) => self._olafft_helper(&mut y, value)
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

    fn _fft_helper(&self, x: &mut [T], h: &mut [T], buffer: &mut [T]) {
        let xlen = x.len();

        let mut xplanner = RealFftPlanner::<T>::new();
        let mut hplanner = RealFftPlanner::<T>::new();
        
        let xfft = xplanner.plan_fft_forward(xlen);
        let hfft = hplanner.plan_fft_forward(xlen);

        let mut xspectrum = xfft.make_output_vec();
        let mut hspectrum = hfft.make_output_vec();

        xfft.process(x, &mut xspectrum).unwrap();
        hfft.process(h, &mut hspectrum).unwrap();

        let mut fft_prod = xspectrum.iter().zip(hspectrum.iter()).map(|(&a, &b)| a * b).collect::<Vec<Complex<T>>>();

        let ifft = xplanner.plan_fft_inverse(xlen);
        let mut ifft_time = ifft.make_output_vec();
        ifft.process(&mut fft_prod, &mut ifft_time).unwrap();

        for (i, value) in buffer.iter_mut().enumerate() {
            *value = ifft_time[i] / T::from(xlen).unwrap(); // return the buffer with len n + m - 1
        }

    }

    fn _olafft_helper(&self, buffer: &mut [T], frame_size: usize) {
        let xlen = self.x.len();
        let hlen = self.h.len();

        let nframes: usize = (xlen as f32 / frame_size as f32).ceil() as usize; // get number of total frames
        let xlen_new = nframes * frame_size; // adjust x len = nframes * frames size
        
        let mut xnewlen = vec![T::zero(); xlen_new];
        xnewlen[..xlen].copy_from_slice(&self.x[..xlen]);

        let len_fft_buffer = frame_size + hlen - 1; // len of fft buffer -> frame size + kernel size - 1
        let ylen = nframes * frame_size + hlen - 1; // len of result vector

        let pow_len = 1 << ((len_fft_buffer as f32).log2() + 1.0) as usize;
        let mut hpad = vec![T::zero(); pow_len];
        hpad[..hlen].copy_from_slice(&self.h[..hlen]);

        let mut xframe = vec![T::zero(); pow_len];

        let mut fft_buffer: Vec<T> = vec![T::zero(); len_fft_buffer]; // fft buffer
        let mut y = vec![T::zero(); ylen]; // result vector
        
        for i in 0..nframes {
            let frame_start = i * frame_size;
            let frame_end = (i + 1) * frame_size;
            xframe[..frame_size].copy_from_slice(&xnewlen[frame_start..frame_end]);
            self._fft_helper(&mut xframe, &mut hpad, &mut fft_buffer);

            for j in 0.. len_fft_buffer {
                y[frame_start + j] += fft_buffer[j]; // rebuild from ola with len nframes * frame size + hsize - 1
            }
        }

        for (i, value) in buffer.iter_mut().enumerate() {
            *value = y[i]; // return the buffer with len n + m - 1
        }

    }

}

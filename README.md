# Convolution Algorithms  

[Def. 1.1]  
$(x*h) [n] = y[n] = \sum_{k=0}^{|h| - 1}h[k] \cdot x[n - k]$

*What algorithms are most time efficient for computing the convolution  
(according to Def. 1.1) of two time-series on a modern machine and in  
which scenarios do they perform best?* [1]  

**Algos**:  

- **Output-side algorithm** $\rightarrow O(|y| \cdot |h|) = O(|x| \cdot |h| + |h|^2)$  
  let $n > |x|, |h| \rightarrow O(n^2)$  
- **Input-side algorithm** $\rightarrow O(|x| \cdot |h|)$  
  let $n > |x|, |h|\rightarrow O(n^2)$
- **FFT convolution** $\rightarrow O(|x| \log{|x|} + |h| \log{|h|})$
- **OLA FFT** $\rightarrow O(k \cdot |n| \log{|n|})$ where k is the number of frames and |n| is the frame size  

**[1]** N. Ghidirimschi, *Convolution Algorithms for interger data types*, University of Groningen, 2021

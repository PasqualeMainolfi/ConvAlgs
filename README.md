# Convolution Algorithms  

[Def. 1.1]  
$(x*h) [n] = y[n] = \sum_{k=0}^{|h| - 1}h[k] \cdot x[n - k]$

*What algorithms are most time efficient for computing the convolution  
(according to Def. 1.1) of two time-series on a modern machine and in  
which scenarios do they perform best?* [1]  

**Algos**:  

- **Output-side algorithm**. The complexity of this algorithm is  
  $O(|y| \cdot |h|) = O(|x| \cdot |h| + |h|^2)$. Let $n > |x|, |h|$, then this can be  
  further expressed as $O(n^2)$  
- **Input-side algorithm**. The complexity of this algorithm is  
  $O(|x| \cdot |h|)$. With $n > |x|, |h|$, the algorithm also belongs to $O(n^2)$
- **Karatsuba (adapted) algorithm**. The complexity of this algorithm is  
  $O(n^{\log_{2}3}) \approx O(n^{1.58})$  
- **FFT convolution** 
- **OLA FFT**

**[1]** N. Ghidirimschi, *Convolution Algorithms for interger data types*, University of Groningen, 2021
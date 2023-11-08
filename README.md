# Decorrelation_Network
Decorrelation Discrete Data Network Algorithm

![mixed_dda](https://user-images.githubusercontent.com/97409887/236999982-bd3f05ca-c3a0-4b37-a7c8-60c40e633a10.png)

Discrete data, often found in medical or survival data, is difficult to work with, especially in the case of correlated data.  To make use of causal inference methods (which use the i.i.d assumption), we de-correlate the discrete data using a novel method that approximates the Expectation-Maximization (EM) Algorithm.  This method shows improved causal graph estimates over baseline methods that assume independent observations.

We use this method on real RNA-seq data and obtain the causal network between different genes in the data.  To determine if our method improves on traditional naive methods of causal inference, we obtain the likelihood of test data through cross-validation and compare between different methods.

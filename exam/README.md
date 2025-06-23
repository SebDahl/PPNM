This folder contains my work for the following task:
Least-squares missing samples recovery (error concealment)
Introduction
Implement the least-squares missing samples recovery method as described in the Book.
Task
Implement a function with the signature
vector recover(vector y)
that takes a corrupted signal y (where the missing elements are replaced by zeros) and returns the reconstructed vector x where the elements that were equal to zero are recovered using the least squares method described in the Book.



My implementation is capable of receiving a signal, finding the missing samples (where all points equal to 0.0 are read as missing) and then recovers the missing samples by minimizing the second derivate using the second order difference matrix.
Based on standard least square statistics I have attempted implemented a calculation of the uncertainty of the recovered samples.

I have added a signal declipping implementation as well, which largely works by the same principle, but finds the clipped samples based on a given threshold and recovers the declipped samples by minimizing the third derivate using the third order difference matrix. This method works fairly well, however it struggles a bit at the boundaries of the samples, but this can be lessened by introducing more samples.


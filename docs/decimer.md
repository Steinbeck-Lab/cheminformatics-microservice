---
outline: deep
---

# ***DECIMER Image Transformer***: Deep Learning for Chemical Image Recognition using Efficient-Net V2 + Transformer

The DECIMER 2.2 [5] (Deep lEarning for Chemical ImagE Recognition) project [1] was launched to address the OCSR problem with the latest computational intelligence methods to provide an automated open-source software solution.

The original implementation of DECIMER[1] using GPU takes a longer training time when we use a bigger dataset of more than 1 million images. To overcome these longer training times, many implement the training script to work on multiple GPUs. However, we tried to step up and implemented our code to use Google's Machine Learning hardware [TPU(Tensor Processing Unit)](https://en.wikipedia.org/wiki/Tensor_Processing_Unit) [2]. You can learn more about the hardware [here](https://en.wikipedia.org/wiki/Tensor_Processing_Unit).

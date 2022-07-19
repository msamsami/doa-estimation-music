# DOA Estimation using the MUSIC algorithm

<p>
<img src="https://www.vocal.com/wp-content/uploads/2019/10/MUSIC-1-768x385.png" width="450" />
<br><br>
</p>

## Introduction
This repository presents a mini project to show how the direction of arrivals (DOA) is estimated using the MUSIC algorithm. MUSIC stands for **MU**ltiple **SI**gnal **C**lassification. It is a subspace-based algorithm used for estimating the DOA of same-frequency narrowband sources arriving at a sensor array. MUSIC outperforms simple methods such as picking peaks of DFT spectra in the presence of noise when the number of components is known in advance because it exploits knowledge of this number to ignore the noise in its final report.

## Code description

In this code, it is assumed that there are *M* sensors with a sampling frequency of 10^7 Hz and *N* independent sources that generate narrowband complex exponential signals with a constant frequency of 10^6 Hz and a standard Gaussian amplitude. For simplicity, it is assumed that the sensors form an array and are positioned in a linear fashion, although more complicated array structures are possible.

## Background
The script in this repository is the completed answer to a course project for the *Linear Algebra* course at *Shiraz University*, which was lectured by [Dr. Alireza Masnadi-Shirazi](https://www.linkedin.com/in/alireza-masnadi-shirazi-aa093561/) back in the Spring 2017 semester. The project guide, a related paper, and a tutorial are also included.

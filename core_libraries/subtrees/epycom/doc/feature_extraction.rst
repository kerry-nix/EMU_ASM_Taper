Feature extraction
======================
This is the core of the whole library. The algorithms for feature extraction are divided into 3 subgroups:

- Univariate
- Bivariate
- Event detection

All the algorithms accept raw or filtered data and provide pandas dataframes as their output.


Univariate feature extraction
*********************************



Bivariate feature extraction
*********************************
Bivariate feature extraction algorithms server for calculating relationships between two signals. They can be used for example to obtain connectivity between different areas of the brain.

- Linear correlation

  The linear correlation (LC) varies in interval <-1,1> and reflects shape similarities between two signals. LC=1 indicates perfect conformity between two signals, LC=-1 indicates opposite signals and LC=0 indicates two different signals. LC is calculated by Pearson’s correlation coefficient as: LCX,Y=[cov(Xt,Yt)/std(Xt)・std(Yt)], where Xt,Yt are the two evaluated signals, cov is the covariance and std is the standard deviation. The linear correlation between two signals can be calculated with a time-lag. Maximum time-lag should not exceed fmax/2. Lagged linear correlation (LLC) for each time-lag k was calculated by Pearson’s correlation coefficient as: LLCX(k),Y(k)=[cov(Xt(k),Yt(k))/std(Xt(k))・std(Yt(k))], where Xt,Yt are the two evaluated signals, cov is the covariance and std is the standard deviation. The maximum value of correlation is stored with its time-lag value.

- Phase consistency

  Phase consistency (PC) varies in interval <0,1> and reflects conformity in phase between two signals, regardless of any phase shift between them. First, phase synchrony (PS) is calculated for multiple steps of time delay between two signals as PS=√[(<cos(ΦXt)>)2+(<sin(ΦYt)>)2], where ΦXt is instantaneous phase of signal X, ΦYt is instantaneous phase of signal Y, <> stands for mean and √ for square root. PC is then calculated as PC = <PS>・(1-std(PS)/0.5), where std is the standard deviation and <> stands for mean.  Instantaneous phase ΦXt is calculated as ΦXt=arctan(xH/xt), where xH is the Hilbert transformation of the time signal xt.

- Phase lag index

  Phase lag index (PLI) varies in interval <0,1> and represents evaluation of statistical interdependencies between time series, which is supposed to be less influenced by the common sources (Stam et al. 2007). PLI calculation is based on the phase synchrony between two signals with constant, nonzero phase lag, which is most likely not caused by volume conduction from a single strong source. Phase lag index is calculated as PLI=|<sign[dΦ(tk)]>|, where sign represents signum function, <> stands for mean and dΦ is a phase difference between two iEEG signals. Maximum time-lag should not exceed fmax/2. The maximum value of PLI is stored with its time-lag value.

- Phase synchrony

  Phase synchrony (PS) varies in interval <0,1> and reflects synchrony in phase between two signals. PS is calculated as PS=√[(<cos(ΦXt)>)2+(<sin(ΦYt)>)2], where ΦXt is instantaneous phase of signal X, ΦYt is instantaneous phase of signal Y, <> stands for mean and √ for square root. Instantaneous phase ΦXt is calculated as ΦXt=arctan(xH/xt), where xH is the Hilbert transformation of the time signal xt.

- Relative entropy

  To evaluate the randomness and spectral richness between two time-series, the Kullback-Leibler divergence, i.e. relative entropy (REN), is calculated. REN is a measure of how entropy of one signal diverges from a second, expected one. The value of REN varies in interval <0,+Inf>. REN=0 indicates the equality of  statistical distributions of two signals, while REN>0 indicates that the two signals are carrying different information. REN is calculated between signals X, Y as REN=sum[pX・log(pX/pY)], where pX is a probability distribution of investigated signal and pY is a probability distributions of expected signal. Because of asymmetrical properties of REN, REN(X, Y) is not equal to REN(Y, X). REN is calculated in two steps for both directions (both distributions from channel pair were used as expected distributions). The maximum value of REN is then considered as the final result, regardless of direction.

- Spectra multiplication

  Frequency spectra of two signals, obtained by Fourier transform, are multiplied and transformed back to the time domain, where the mean and std is calculated.


Event detection
*********************************
This subsection provides algorithms for detection of events occurring in the signal. All algorithms provide event position or event start/stop and some of them provide additional features of detected events. Currently the library contains algorithms for detecting interictal epileptiform discharges (IEDs),i.e. epileptic spikes, and a number of algorithms for detection of high frequency oscillations (HFOs).

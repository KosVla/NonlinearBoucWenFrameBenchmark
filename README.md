# Two-story frame with Bouc-Wen hysteretic links as a multi-degree of freedom nonlinear response simulator

## Code repository
The software implementation of a two-story frame with Bouc-Wen hysteretic links is provided in this repository, as part of a multi-degree of freedom nonlinear response simulator benchmark case study proposed in the ** 5th Edition of the Workshop on Nonlinear System Identification Benchmarks** (April, 2021, [Link](https://sites.google.com/view/nonlinear-benchmark/benchmarks)). Together with the software implementation, the benchmark simulator has been used to create five standardized datasets for identification applications. Those are provided through the zenodo platform hereLink for reference purposes.

The implementation of the software aims to be utilized as a benchmark problem to validate methods and tools in structural health monitoring, model reduction, or identification applications. The proposed benchmark is provided as a framework simulator and not a single function, thus offering full flexibility for the user to modify and evaluate the shear frame based on customized needs and requirements of the underlying problem. For this reason, the frame is excited using a parametrized ground motion excitation, created through a stochastic signal, the Bouc-Wen model is parametrized, and the software offers extension possibilities such as multi-story frame assembly, deterioration or degradation phenomena, and localized damage representation. All these possibilities are documented in the description provided. 

The main branch contains the description of the case study, along with a welcome *ReadMe*.
The standardized MATLAB implementation is provided in the *Version0.0_Matlab* branch, along with a dedicated *ReadMe* file documenting the coding scripts and the subdirectories of the repository. In this branch the files used to simulate the standardized datasets are also provided, along with an additional *ReadMe* documenting the configuration files.

An example python implementation of the framework is also provided in the *Version0.0_Python* branch, along with a dedicated *ReadMe* file.

The additional branches provide updated or modified versions of the software. Each branch contains a dedicated *ReadMe* file documenting the adjustments compared to the standardized version.

## Motivation
A diverse variety of engineering and dynamic systems, ranging from control applications and solid mechanics to biology and economics, feature hysteretic phenomena. This commonly encountered nonlinear behavior can be captured and described via diverse numerical models, with the Bouc-Wen representation comprising a common choice within the nonlinear dynamics and vibration engineering community. In this benchmark, the Bouc-Wen model is employed to characterize the response of the nodal connections of a two-story frame structure. The resulting shear frame with hysteretic links is proposed as a multi-degree of freedom nonlinear response simulator.

This case study can be seen as an extension to the single degree of freedom 'Hysteretic Benchmark with a Dynamic Nonlinearity' problem, which is already featured in the nonlinear benchmark catalogue [here](https://sites.google.com/view/nonlinear-benchmark/). Our simulator employs a similar parameterized representation of the Bouc-Wen model for each nonlinear link and builds upon it to also include strength deterioration and stiffness degradation effects in a structure with increased dimensionality. As a result, the featured parametric shear frame serves as a multi-degree of freedom nonlinear response simulator, able to model a wide range of nonlinear effects through the parametrized Bouc-Wen couplings.

The provided simulator can be utilized as a benchmark problem to validate methods and tools in structural health monitoring, model reduction, or identification applications. It has already been exploited in a reduced order modelling context to validate the performance of parametric, physics-based ROMs for nonlinear, dynamical systems in [2,3,4]. 
Compared to the existing Bouc-Wen oscillator benchmark featured in the nonlinear benchmark website, described in detail in [5], the proposed multi-degree of freedom simulator allows for increased complexity studies due to the higher dimensionality of the system and the potential for multi-parametric numerical examples. In addition, the proposed benchmark is provided as a framework simulator and not a single function, thus offering full flexibility for the user to modify and evaluate the shear frame based on customized needs and requirements of the underlying problem.

## References
[1] J.P. Noel, M. Schoukens, and K. Tiels, “Nonlinear benchmarks,” https://sites.google.com/view/nonlinear-benchmark/.

[2] K. Vlachas, K. Tatsis, K. Agathos, A. R. Brink, and E. Chatzi, “A physics-based, local POD basis approach formulti-parametric reduced order models,” in International Conference on Noise and Vibration Engineering (ISMA2020) in conjunction with the 8th International Conference on Uncertainty in Structural Dynamics (USD 2020).ETH Zurich, Environmental and Geomatic Engineering, 2020.

[3] K. Vlachas, K. Tatsis, K. Agathos, A. Brink, and E. Chatzi, “A local basis approximation approach for nonlinearparametric model order reduction,”Journal of Sound and Vibration, vol. 502, p. 116055, 2021.

[4]  T. Simpson,  N. Dervilis,  and E. Chatzi,  “A machine learning approach to model order reduction of nonlinear systems via autoencoder and LSTM networks,” Journal of Engineering Mechanics, (Forthcoming)

[5] J. P. No ̈el and M. Schoukens, “Hysteretic benchmark with a dynamic nonlinearity,” in Workshop on nonlinearsystem identification benchmarks, 2016, pp. 7–14.

[6] R. Bouc, “A mathematical model for hysteresis (Mod`ele math`ematique d’hyst`er`esis: application aux syst`emes`aun degr`e de libert`e),” Acustica (in French), vol. 24, pp. 16–25, 1971.

[7] R. Bouc, “Forced vibrations of mechanical systems with hysteresis,”Proc. of the Fourth Conference on Nonlinear Oscillations, Prague, 1967.

[8] Y.-K. Wen, “Method for random vibration of hysteretic systems,” Journal of the engineering mechanics division,vol. 102, no. 2, pp. 249–263, 1976.

[9] F. Ikhouane and J. Rodellar, Systems with hysteresis: analysis, identification and control using the Bouc-Wen model.    John Wiley & Sons, 2007.

[10] M. Ismail, F. Ikhouane, and J. Rodellar, “The hysteresis bouc-wen model, a survey,” Archives of Computational Methods in Engineering, vol. 16, no. 2, pp. 161–188, 2009.

[11] T. T. Baber and Y.-K. Wen, “Random vibration of hysteretic, degrading systems,” Journal of the Engineering Mechanics Division, vol. 107, no. 6, pp. 1069–1087, 1981.

## Description
The attached description document provides a detailed description of the proposed multi-degree of freedom nonlinear response simulator. User guidelines for proper input and accurate simulations are also provided, along with a detailed overview of the capabilities of the provided software. This serves as a short overview of the user to treat the proposed benchmark as a fully adjustable case study. The challenges associated with performing nonlinear system identification in the featured benchmark are discussed and standardized system identification tasks related to the benchmark along with the respective provided datasets are proposed. 

## Features

* Multi degree of freedom nonlinear response simulator
* Hysteretic behavior of connections through Bouc-Wen links
* Degradation and deterioration phenomena
* Stochastic ground motion excitation
* User-input excitation (signal) possible
* Parametric dynamic response and parametrized Bouc-Wen links
* Multi-story extension possible through automatic input file creation function  
* Simulation of localized phenomena
* Standardized datasets for tasks related to system identification applications, reduced-order or surrogate modelling applications. 

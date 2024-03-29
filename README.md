# Two-story frame with Bouc-Wen hysteretic links as a multi-degree of freedom nonlinear response simulator

## Code repository
The software implementation of a two-story frame with Bouc-Wen hysteretic links is provided in this repository, as part of a multi-degree of freedom nonlinear response simulator benchmark case study proposed in the ** 5th Edition of the Workshop on Nonlinear System Identification Benchmarks** (April, 2021, [Link](https://sites.google.com/view/nonlinear-benchmark/benchmarks)). Together with the software implementation, the benchmark simulator has been used to create five standardized datasets for identification applications. Those are provided through the zenodo platform [here](https://doi.org/10.5281/zenodo.4742248) for reference purposes.

The implementation of the software aims to be utilized as a benchmark problem to validate methods and tools in structural health monitoring, model reduction, or identification applications. The proposed benchmark is provided as a framework simulator and not a single function, thus offering full flexibility for the user to modify and evaluate the shear frame based on customized needs and requirements of the underlying problem. For this reason, the frame is excited using a parametrized ground motion excitation, created through a stochastic signal, the Bouc-Wen model is parametrized, and the software offers extension possibilities such as multi-story frame assembly, deterioration or degradation phenomena, and localized damage representation. All these possibilities are documented in the description provided. 

## Dataset Repository 

The benchmark simulator has been used to derive five standardized simulation datasets as part of the conference submission to be employed as a comparison reference for identification applications. Each dataset refers to a separate scenario, representing damage phenomena, parameter estimation or calibration case studies. The detailed setup of the datasets can be found on the description file. The datasets are provided through the zenodo platform [here](https://doi.org/10.5281/zenodo.4742248).

## Repository structure

The **main** branch contains the description of the case study, along with a welcome *ReadMe* file.

The standardized MATLAB implementation is provided in the **Version0.0_Matlab** branch, along with a dedicated *ReadMe* file documenting the coding scripts and the subdirectories of the repository. In this branch the files used to simulate the standardized datasets are also provided, along with an additional *ReadMe* documenting the configuration files.

An example python implementation of the framework is also provided in the **Version0.0_Python** branch, along with a dedicated *ReadMe* file.

Any additional branches provide updated or modified versions of the software. Each branch contains a dedicated *ReadMe* file documenting the adjustments compared to the standardized version. 

For example, the  **Version_1.0_Matlab** branch contains a modified frame with hysteretic links assembled only on the rotational degrees of freedom and bending-to-shear coupling to represent a shear frame under earthquake excitation in a more realistic way. All corresponding modifications are documented on the respective *ReadMe* files of the branch.
Modified Features compared to Version 0.0 include:
* Bouc-Wen links are activated only on the rotational degrees of freedom
* Bending and shear degrees of freedom are now coupled
* The option to assemble the Bouc-Wen bw_k coefficient from the EI coefficient of the respective beam element


## Features

* Multi degree of freedom nonlinear response simulator
* Hysteretic behavior of connections through Bouc-Wen links
* Degradation and deterioration phenomena
* Stochastic ground motion excitation
* User-input excitation (signal) possible
* Parametric dynamic response and parametrized Bouc-Wen links
* Multi-story extension possible through automatic input file creation function  
* Simulation of localized phenomena
* Standardized datasets for tasks related to system identification applications, reduced-order or surrogate modelling applications in [here](https://doi.org/10.5281/zenodo.4742248)

## Motivation
A diverse variety of engineering and dynamic systems, ranging from control applications and solid mechanics to biology and economics, feature hysteretic phenomena. This commonly encountered nonlinear behavior can be captured and described via diverse numerical models, with the Bouc-Wen representation comprising a common choice within the nonlinear dynamics and vibration engineering community. In this benchmark, the Bouc-Wen model is employed to characterize the response of the nodal connections of a two-story frame structure. The resulting shear frame with hysteretic links is proposed as a multi-degree of freedom nonlinear response simulator.

This case study can be seen as an extension to the single degree of freedom 'Hysteretic Benchmark with a Dynamic Nonlinearity' problem, which is already featured in the nonlinear benchmark catalogue [here](https://sites.google.com/view/nonlinear-benchmark/). Our simulator employs a similar parameterized representation of the Bouc-Wen model for each nonlinear link and builds upon it to also include strength deterioration and stiffness degradation effects in a structure with increased dimensionality. As a result, the featured parametric shear frame serves as a multi-degree of freedom nonlinear response simulator, able to model a wide range of nonlinear effects through the parametrized Bouc-Wen couplings.

The provided simulator can be utilized as a benchmark problem to validate methods and tools in structural health monitoring, model reduction, or identification applications. It has already been exploited in a reduced order modelling context to validate the performance of parametric, physics-based ROMs for nonlinear, dynamical systems in [2,3,4]. 
Compared to the existing Bouc-Wen oscillator benchmark featured in the nonlinear benchmark website, described in detail in [5], the proposed multi-degree of freedom simulator allows for increased complexity studies due to the higher dimensionality of the system and the potential for multi-parametric numerical examples. In addition, the proposed benchmark is provided as a framework simulator and not a single function, thus offering full flexibility for the user to modify and evaluate the shear frame based on customized needs and requirements of the underlying problem.

## References
[1] J.P. Noël, M. Schoukens, and K. Tiels, “Nonlinear benchmarks,” https://sites.google.com/view/nonlinear-benchmark/.

[2] K. Vlachas, K. Tatsis, K. Agathos, A. R. Brink, and E. Chatzi, “A physics-based, local POD basis approach for multi-parametric reduced order models,” in International Conference on Noise and Vibration Engineering (ISMA2020) in conjunction with the 8th International Conference on Uncertainty in Structural Dynamics (USD 2020).

[3] K. Vlachas, K. Tatsis, K. Agathos, A. Brink, and E. Chatzi, “A local basis approximation approach for nonlinear parametric model order reduction,” Journal of Sound and Vibration, vol. 502, p. 116055, 2021.

[4]  T. Simpson,  N. Dervilis,  and E. Chatzi,  “A machine learning approach to model order reduction of nonlinear systems via autoencoder and LSTM networks,” Journal of Engineering Mechanics, (Forthcoming)

[5] J. P. Noël and M. Schoukens, “Hysteretic benchmark with a dynamic nonlinearity,” in Workshop on nonlinear system identification benchmarks, 2016, pp. 7–14.

[6] R. Bouc, “A mathematical model for hysteresis,” Acustica (in French), vol. 24, pp. 16–25, 1971.

[7] R. Bouc, “Forced vibrations of mechanical systems with hysteresis,” Proc. of the Fourth Conference on Nonlinear Oscillations, Prague, 1967.

[8] Y.-K. Wen, “Method for random vibration of hysteretic systems,” Journal of the engineering mechanics division,vol. 102, no. 2, pp. 249–263, 1976.

[9] F. Ikhouane and J. Rodellar, Systems with hysteresis: analysis, identification and control using the Bouc-Wen model.    John Wiley & Sons, 2007.

[10] M. Ismail, F. Ikhouane, and J. Rodellar, “The hysteresis bouc-wen model, a survey,” Archives of Computational Methods in Engineering, vol. 16, no. 2, pp. 161–188, 2009.

[11] T. T. Baber and Y.-K. Wen, “Random vibration of hysteretic, degrading systems,” Journal of the Engineering Mechanics Division, vol. 107, no. 6, pp. 1069–1087, 1981.

## Description
The attached description document provides a detailed description of the proposed multi-degree of freedom nonlinear response simulator. User guidelines for proper input and accurate simulations are also provided, along with a detailed overview of the capabilities of the provided software. This serves as a short overview of the user to treat the proposed benchmark as a fully adjustable case study. The challenges associated with performing nonlinear system identification in the featured benchmark are discussed and standardized system identification tasks related to the benchmark along with the respective provided datasets are proposed. 

## Works where the benchmark is featured

Previous versions of the benchmark simulator are featured in the following works:

Simpson, Thomas, Nikolaos Dervilis, and Eleni Chatzi. "Machine Learning Approach to Model Order Reduction of Nonlinear Systems via Autoencoder and LSTM Networks." Journal of Engineering Mechanics 147, no. 10 (2021): 04021061.

Vlachas, Konstantinos, Konstantinos Tatsis, Konstantinos Agathos, Adam R. Brink, and Eleni Chatzi. "A local basis approximation approach for nonlinear parametric model order reduction." Journal of Sound and Vibration 502 (2021): 116055.

Vlachas, Konstantinos, Konstantinos Tatsis, Konstantinos Agathos, Adam R. Brink, and Eleni Chatzi. "A physics-based, local POD basis approach for multi-parametric reduced order models." In International Conference on Noise and Vibration Engineering (ISMA 2020) in conjunction with the 8th International Conference on Uncertainty in Structural Dynamics (USD 2020). ETH Zurich, Environmental and Geomatic Engineering, 2020.

## Cite

If you use the algorithmic framework presented here you are kindly requested to cite the following work(s):

```
@inproceedings{vlachas2021two,
  title={Two-story frame with Bouc-Wen hysteretic links as a multi-degree of freedom nonlinear response simulator},
  author={Vlachas, Konstantinos and Agathos, Konstantinos and Tatsis, Konstantinos E and Brink, Adam R and Chatzi, Eleni},
  booktitle={5th Workshop on Nonlinear System Identification Benchmarks (2021)},
  pages={6},
  year={2021}
}
```
}

Further sources which can be of interest include:
```
@article{vlachas2021local,
  title={A local basis approximation approach for nonlinear parametric model order reduction},
  author={Vlachas, Konstantinos and Tatsis, Konstantinos and Agathos, Konstantinos and Brink, Adam R and Chatzi, Eleni},
  journal={Journal of Sound and Vibration},
  volume={502},
  pages={116055},
  year={2021},
  publisher={Elsevier}
}
```

```
@article{simpson2021machine,
  title={Machine learning approach to model order reduction of nonlinear systems via autoencoder and LSTM networks},
  author={Simpson, Thomas and Dervilis, Nikolaos and Chatzi, Eleni},
  journal={Journal of Engineering Mechanics},
  volume={147},
  number={10},
  pages={04021061},
  year={2021},
  publisher={American Society of Civil Engineers}
}
```

```
@inproceedings{vlachas2022coupling,
  title={On the coupling of reduced order modeling with substructuring of structural systems with component nonlinearities},
  author={Vlachas, Konstantinos and Tatsis, Konstantinos and Agathos, Konstantinos and Brink, Adam R and Quinn, Dane and Chatzi, Eleni},
  booktitle={Dynamic Substructures, Volume 4: Proceedings of the 39th IMAC, A Conference and Exposition on Structural Dynamics 2021},
  pages={35--43},
  year={2022},
  organization={Springer}
}
```

# Shear frame with Bouc-Wen hysteretic links as a multi-degree of freedom nonlinear response simulator

The numerical implementation of a shear structural frame with Bouc-Wen hysteretic links is featured here as a benchmark problem to validate methods and tools in structural health monitoring, model reduction, or structural system identification applications. It has been originally presented in the **5th Edition of the Workshop on Nonlinear System Identification Benchmarks** (April, 2021, [Link](https://sites.google.com/view/nonlinear-benchmark/benchmarks)). 

The proposed benchmark is provided as a framework simulator and not a single function, thus offering full flexibility for the user to modify and evaluate the shear frame based on customized needs and requirements of the underlying problem. For this reason, the frame is excited using a parametrized ground motion excitation, created through a stochastic signal, the Bouc-Wen model is parametrized, and the software offers extension possibilities such as multi-story frame assembly, deterioration or degradation phenomena, and localized damage representation. All these possibilities are documented in the description provided. 

## Features

* Multi degree of freedom nonlinear response simulator
* Hysteretic behavior of connections through Bouc-Wen links
* Degradation and deterioration phenomena
* Stochastic ground motion excitation
* User-input excitation (signal) possible
* Parametric dynamic response and parametrized Bouc-Wen links
* Multi-story extension possible through automatic input file creation function  
* Simulation of localized phenomena
* Standardized datasets for benchmarking in [here](https://doi.org/10.5281/zenodo.4742248)

## Code repository
The software implementation is provided in this repository in both `MATLAB` and `Python`.

## Dataset for structural system identification tasks

The benchmark simulator has been used to derive five standardized simulation datasets to be employed as a comparison reference for system identification applications. Each dataset refers to a separate scenario, representing damage phenomena, parameter estimation or calibration case studies. The detailed setup of the datasets can be found on the description file. The datasets are provided through the zenodo platform [here](https://doi.org/10.5281/zenodo.4742248).

## Repository structure

The **main** branch contains the description of the case study, along with a welcome *ReadMe* file. The standardized MATLAB implementation is provided in the **MatlabImplementation** directory, along with a dedicated *ReadMe* file. A python implementation of the framework is also provided in the **PythonImplementation** directory, along with a dedicated *ReadMe* file.

## Works where the benchmark is featured

The benchmark simulator is featured in the following works:

Simpson, Thomas, et al. "VpROM: a novel variational autoencoder-boosted reduced order model for the treatment of parametric dependencies in nonlinear systems." Scientific Reports 14.1 (2024): 6091.[DOI](https://doi.org/10.1038/s41598-024-56118-x)

Vlachas, Konstantinos, et al. "Parametric reduced-order modeling for component-oriented treatment and localized nonlinear feature inclusion." Nonlinear Dynamics 112.5 (2024): 3399-3420.[DOI](https://doi.org/10.1007/s11071-023-09213-z)

Vlachas, Konstantinos, et al. "Reduced Order Modeling conditioned on monitored features for response and error bounds estimation in engineered systems." Mechanical Systems and Signal Processing 226 (2025): 112261.[DOI](https://doi.org/10.1016/j.ymssp.2024.112261)

Kamariotis, Antonios, et al. "On the consistent classification and treatment of uncertainties in structural health monitoring applications." ASCE-ASME Journal of Risk and Uncertainty in Engineering Systems, Part B: Mechanical Engineering 11.1 (2025): 011108.[DOI](https://doi.org/10.1115/1.4067140)



## Cite

If you use the algorithmic framework presented here you are kindly requested to cite the following work(s):

```bibtex
@article{Vlachas2021,
  title={A local basis approximation approach for nonlinear parametric model order reduction},
  author={Vlachas, Konstantinos and Tatsis, Konstantinos and Agathos, Konstantinos and Brink, Adam R and Chatzi, Eleni},
  journal={Journal of Sound and Vibration},
  volume={502},
  pages={116055},
  year={2021},
  publisher={Elsevier},
  doi={10.1016/j.jsv.2021.116055}
}
```

```bibtex
@inproceedings{vlachas2021two,
  title={Two-story frame with Bouc-Wen hysteretic links as a multi-degree of freedom nonlinear response simulator},
  author={Vlachas, Konstantinos and Agathos, Konstantinos and Tatsis, Konstantinos E and Brink, Adam R and Chatzi, Eleni},
  booktitle={5th Workshop on Nonlinear System Identification Benchmarks (2021)},
  pages={6},
  year={2021}
}
```


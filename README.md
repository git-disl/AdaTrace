# AdaTrace

This repository contains a snapshot of the code for AdaTrace: utility-aware, differentially private, and attack-resilient location trace synthesizer. Our paper describing the AdaTrace system was published in ACM CCS 2018:

```
Utility-Aware Synthesis of Differentially Private and Attack-Resilient Location Traces. 
Gursoy, M. E., Liu, L., Truex, S., Yu, L., & Wei, W. (2018, October). 
In Proceedings of the 2018 ACM SIGSAC Conference on Computer and Communications Security (pp. 196-211). ACM.
```

## Contents 

Source code can be found under the `src` folder. 

For ease of use, we have prepared an additional readme file titled `AdaTrace-Readme.docx` which contains detailed screenshots and explanations of the trajectory synthesis and experimentation process. 

We also provide a sample dataset: `brinkhoff.dat`. 

## Installation and Dependencies

AdaTrace source code is written in Java. There are two dependencies on external libraries: `commons-math3-3.4.1.jar` and `kd.jar`. Both are provided in this repository. Please add them to your project build path before trying to build AdaTrace. 

## Project Status

We are continuing the development of AdaTrace and there is ongoing work in our lab regarding privacy-preserving trajectory analytics. Meanwhile, we provide this code for quick dissemination of research artifacts. 

The code is provided as is, without warranty or support. If you use our code, please cite:

```
@inproceedings{gursoy2018utility,
  title={Utility-Aware Synthesis of Differentially Private and Attack-Resilient Location Traces},
  author={Gursoy, Mehmet Emre and Liu, Ling and Truex, Stacey and Yu, Lei and Wei, Wenqi},
  booktitle={Proceedings of the 2018 ACM SIGSAC Conference on Computer and Communications Security},
  pages={196--211},
  year={2018},
  organization={ACM}
}
```

```
@article{gursoy2018differentially,
  title={Differentially private and utility preserving publication of trajectory data},
  author={Gursoy, Mehmet Emre and Liu, Ling and Truex, Stacey and Yu, Lei},
  journal={IEEE Transactions on Mobile Computing},
  year={2018},
  publisher={IEEE}
}
```

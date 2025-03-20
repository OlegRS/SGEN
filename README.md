# SGEN (Stochastic Gene Expression in Neurons)

SGEN is a high-performance library for computing distributions of **mRNAs and proteins in arbitrary neurons**. The core of the library is written in **C++** and wrapped in **Python** using `pybind11`. 

SGEN is primarily designed to be used through the `SGEN_Py` **Python module** but can also be used directly in C++.

## Installation

### Local Installation
To install the library on your system, run:
```bash
git clone https://github.com/OlegRS/SGEN && cd SGEN && ./setup.sh
```

### Google Colab Installation
To use the library in **Google Colab**, insert the following commands into a code cell:
```bash
!git clone https://github.com/OlegRS/SGEN
%cd SGEN
!apt install -qq xvfb libgl1-mesa-glx
!./setup.sh
```

## Usage
### Quick Start
For an interactive **tutorial** on using SGEN, check out the **Google Colab Quick Start**:
[Open Tutorial in Colab](https://colab.research.google.com/drive/1MTaI-sJZjxmD___74ecXpNLfsrD2qbRW?usp=sharing)

### Examples
- **Python usage examples**: [`SGEN/python_examples`](https://github.com/OlegRS/SGEN/tree/main/python_examples)
- **C++ usage examples**: [`SGEN/cpp_examples`](https://github.com/OlegRS/SGEN/tree/main/cpp_examples)

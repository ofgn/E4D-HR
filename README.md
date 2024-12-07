
# E4D-HR

E4D-HR is a fork of [E4D](https://www.pnnl.gov/projects/e4d), a 3D geophysical modelling and inversion tool originally developed by PNNL. This project aims to continue the development of E4D, introducing modifications tailored towards mineral exploration use cases. E4D is designed for subsurface imaging, using direct current resistivity and induced polarisation methods for geophysical modelling and inversion.

---

## Installation

### System Requirements

E4D-HR is compatible with the following Linux distributions and Windows systems that support WSL2 (Windows Subsystem for Linux):

- **Linux**: Ubuntu 22.04, Ubuntu 24.04, Debian 11, Debian 12
- **Windows (WSL2)**: Windows 10 (Version 1903 or later), Windows 11

### Dependencies

Ensure the following software dependencies are installed:

- [Intel® Fortran Essentials](https://www.intel.com/content/www/us/en/developer/tools/oneapi/hpc-toolkit-download.html?packages=fortran-essentials&fortran-essentials-os=linux&fortran-essentials-lin=offline)

  Subset of Intel® oneAPI HPC toolkit
  
- [Intel® C++ Essentials](https://www.intel.com/content/www/us/en/developer/tools/oneapi/base-toolkit-download.html?packages=cpp-essentials&cpp-essentials-os=linux&cpp-essentials-lin=offline)

  Subset of Intel® oneAPI Base toolkit
  
- [build-essential](https://packages.ubuntu.com/oracular/build-essential)
  
  Install through APT
  ```bash
  sudo apt-get update
  sudo apt-get install build-essential
  ```

### Download

1. Visit the [E4D-HR GitHub repository](https://github.com/ofgn/E4D-HR).
   
2. Click the **Code** button on the repository page.
   
3. Choose one of the following methods:
   
   - Download the repository as a ZIP file and extract it to your desired directory.
     
   - Clone the repository with Git:
     
     ```bash
     git clone https://github.com/ofgn/E4D-HR.git
     ```
   

### Automatic Installation

To automatically install E4D-HR, use the provided `install.py` script. This script handles the setup of dependencies and the compilation of E4D-HR.

1. Ensure the Intel oneAPI environment variables environment variables are configured.

   [Use the setvars and oneapi-vars Scripts with Linux](https://www.intel.com/content/www/us/en/docs/oneapi/programming-guide/2025-0/use-the-setvars-and-oneapi-vars-scripts-with-linux.html#HOW-TO-RUN)

   ```bash
     source <oneapi-dir>/setvars.sh
     ```
   
3. Navigate to the E4D-HR directory containing the `install.py` script.
   
4. Run the script:
   ```python3 install.py```

The script will build E4D-HR and third party libraries, prompting as needed.

### Manual Installation

See the E4D-HR user guide for information regarding manual installation.

---

## Contributing

E4D-HR is maintained by [ExploreGeo](https://www.exploregeo.com.au/), a consulting group specialising in geophysical services in the minerals exploration industry.

<p align="center">
  <a href="https://www.exploregeo.com.au" target="_blank">
    <img src="https://www.exploregeo.com.au/images/exploregeo/EG_sign.png" alt="EG_sign" />
  </a>
</p>

Contributions to E4D-HR are welcomed! To contribute, please fork the repository, make your changes, and submit a pull request. For any questions or suggestions, feel free to open an issue on GitHub.


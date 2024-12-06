
# E4D-HR

E4D-HR is a fork of [E4D](https://www.pnnl.gov/projects/e4d), a 3D geophysical modelling and inversion tool originally developed by PNNL. This project continues the development of E4D, introducing modifications tailored towards mineral exploration use cases. E4D is designed for subsurface imaging, using direct current resistivity and induced polarisation methods for advanced geophysical modelling and inversion.

---

## Contributing

E4D-HR is maintained by [ExploreGeo](https://www.exploregeo.com.au/), a consulting group specialising in geophysical services tailored to the minerals exploration industry.

Contributions to E4D-HR are welcomed! To contribute, please fork the repository, make your changes, and submit a pull request. For any questions or suggestions, feel free to open an issue on GitHub.

---

## Installation

### 1. System Requirements

E4D-HR is compatible with the following Linux distributions and Windows systems that support WSL2 (Windows Subsystem for Linux):

- **Linux**: Ubuntu 20.04, 22.04, and 24.04
- **Windows (WSL2)**: Windows 10 (Version 1903 or later), Windows 11

### 2. Required Software

Ensure the following software dependencies are installed:

- **Intel oneAPI Base Toolkit**: Version 2024.2.0 or later
- **Intel HPC Toolkit**: Version 2024.2.0 or later
- **build-essential** (for Ubuntu and Debian):
  ```bash
  sudo apt-get update && sudo apt-get install build-essential
  ```

### 3. Downloading E4D-HR

1. Visit the [E4D-HR GitHub repository](https://github.com/ofgn/E4D-HR).
2. Click the **Code** button on the repository page.
3. Choose one of the following methods:
   - Clone the repository with Git:
     ```bash
     git clone https://github.com/ofgn/E4D-HR.git
     ```
   - Download the repository as a ZIP file and extract it to your desired directory.

### 4. Automatic Installation

To automatically install E4D-HR, use the provided `install.py` script. This script handles the setup of dependencies and the compilation of E4D-HR.

1. Ensure the Intel MKL environment variables are configured.
2. Navigate to the directory containing the `install.py` script.
3. Run the script:
   ```python3 install.py```

The script will automatically handle dependency installation and E4D compilation, prompting as needed.

### 5. Manual Installation

See the E4D-HR user guide for information regarding manual installation.


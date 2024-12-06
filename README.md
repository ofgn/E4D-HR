# E4D-HR

E4D-HR is a fork of **E4D**, a 3D geophysical modelling and inversion tool originally developed by PNNL. E4D-HR includes modifications targeting mineral exploration use cases, quality-of-life improvements, and critical bug fixes. **E4D** is a 3D geophysical modeling and inversion code designed for subsurface imaging using Direct Current Resistivity (DCR) and Induced Polarisation (IP) methods.

---

## Key Features of E4D-HR

### 1. **VTK Visualisation**
   - Enhanced visualisation capabilities using VTK for efficient and intuitive inspection of 3D geophysical data.

### 2. **Performance**
   - Faster and more efficient finite element solver for improved computational performance.

### 3. **Bug Fixes**
   - Addressed known issues from the original E4D to improve stability and reliability.

### 4. **Improved Starting Models**
   - Added support for initialising inversion processes with **mean** or **median** models for greater flexibility.

---

## Potential Future Additions

- **Wavelet Compression for the Jacobian Matrix**  
  Efficient compression techniques to optimise memory usage and computation during inversion.

- **Revamped Input File Formats**  
  Simplified and modernised formats for easier setup and greater flexibility.

- **No Third Party Dependencies**  
  Ideally replace all third party libraries with native Fortran code.

---

## Maintenance

E4D-HR is maintained by **ExploreGeo**, a consulting group specialising in geophysical services tailored to the minerals exploration industry. 

For more information about ExploreGeo, visit their [website](https://www.exploregeo.com.au/).

---

## System Requirements

E4D-HR supports the following platforms:

- **Ubuntu**: 20.04 (Focal Fossa), 22.04 (Jammy Jellyfish), 24.04 (Noble Numbat)
- **WSL2**: Windows 10 (Version 1903 or later), Windows 11

---


## Installation

For installation instructions, refer to the [documentation](#). Ensure all required dependencies and software prerequisites are installed prior to setup.

---

## Contributing

Contributions to E4D-HR are welcomed!

---

## Licence

E4D-HR is released under the [E4D licence](#). Ensure you review the terms and conditions before using the software.

---

For further details, support, and technical guidance, visit the [E4D-HR Documentation](#).

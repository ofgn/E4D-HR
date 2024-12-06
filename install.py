import os
import subprocess
import shutil

# Configuration
CONFIG = {
    "software": {
        "petsc_archive": "petsc-3.21.4.tar.gz",
        "petsc_dir_name": "petsc-3.21.4",
        "triangle_archive": "triangle-1.6.tar.gz",
        "triangle_dir_name": "triangle-1.6",
        "tetgen_archive": "tetgen-1.6.0.tar.gz",
        "tetgen_dir_name": "tetgen-1.6.0",
    },
    "compiler": {
        "CC": "icx",
        "CXX": "icpx",
        "FC": "ifx",
        "MPIFC": "mpiifort -fc=ifx",
        "MPICC": "mpiicc -cc=icx",
        "MPICXX": "mpiicpc -cxx=icpx",
        "FOPTFLAGS": "-O2 -march=native -mtune=native",
        "COPTFLAGS": "-O2 -march=native -mtune=native",
        "CXXOPTFLAGS": "-O2 -march=native -mtune=native",
    },
    "directories": {
        "e4d_dir": os.getcwd(),
        "lib_dir": os.path.join(os.getcwd(), "lib"),
        "bin_dir": os.path.join(os.getcwd(), "bin"),
    },
}

# Helper Functions
def is_installed(path):
    """Check if a directory or file exists."""
    return os.path.exists(path)

def run_command(cmd, cwd=None, env=None):
    """Run a shell command with an optional environment."""
    subprocess.run(cmd, shell=True, cwd=cwd, env=env, check=True)

def check_archive_exists(archive_path):
    """Check if the archive file exists in the lib directory."""
    if not os.path.isfile(archive_path):
        print(f"Error: Archive '{archive_path}' not found.")
        print("Ensure the archive is located in the 'lib' directory.")
        exit(1)

def prompt_reinstall():
    """Ask if the user wants to reinstall each component."""
    components = {
        "PETSc": os.path.join(CONFIG["directories"]["lib_dir"], CONFIG["software"]["petsc_dir_name"]),
        "Triangle": os.path.join(CONFIG["directories"]["lib_dir"], CONFIG["software"]["triangle_dir_name"]),
        "TetGen": os.path.join(CONFIG["directories"]["lib_dir"], CONFIG["software"]["tetgen_dir_name"]),
        "E4D-HR": os.path.join(CONFIG["directories"]["bin_dir"], "e4d"),
    }
    reinstall = {}
    for name, path in components.items():
        if is_installed(path):
            response = input(f"{name} is already installed. Reinstall? [y/N]: ").strip().lower()
            reinstall[name] = response == "y"
        else:
            reinstall[name] = True  # Install if not already installed
    return reinstall

def install_petsc():
    """Install PETSc."""
    petsc_archive = os.path.join(CONFIG["directories"]["lib_dir"], CONFIG["software"]["petsc_archive"])
    petsc_dir = os.path.join(CONFIG["directories"]["lib_dir"], CONFIG["software"]["petsc_dir_name"])
    check_archive_exists(petsc_archive)
    shutil.rmtree(petsc_dir, ignore_errors=True)
    print(f"Installing PETSc from {petsc_archive}...")
    shutil.unpack_archive(petsc_archive, CONFIG["directories"]["lib_dir"])
    run_command(
        f"./configure --with-debugging=0 "
        f"--with-fc='{CONFIG['compiler']['MPIFC']}' --with-cc='{CONFIG['compiler']['MPICC']}' "
        f"--with-cxx='{CONFIG['compiler']['MPICXX']}' --with-blaslapack-dir='{os.getenv('MKLROOT')}' "
        f"FOPTFLAGS='{CONFIG['compiler']['FOPTFLAGS']}' COPTFLAGS='{CONFIG['compiler']['COPTFLAGS']}' "
        f"CXXOPTFLAGS='{CONFIG['compiler']['CXXOPTFLAGS']}'",
        cwd=petsc_dir,
    )
    run_command(f"make PETSC_DIR={petsc_dir} PETSC_ARCH=arch-linux-c-opt all", cwd=petsc_dir)
    run_command(f"make PETSC_DIR={petsc_dir} PETSC_ARCH=arch-linux-c-opt check", cwd=petsc_dir)

def install_triangle():
    """Install Triangle."""
    triangle_archive = os.path.join(CONFIG["directories"]["lib_dir"], CONFIG["software"]["triangle_archive"])
    triangle_dir = os.path.join(CONFIG["directories"]["lib_dir"], CONFIG["software"]["triangle_dir_name"])
    check_archive_exists(triangle_archive)
    shutil.rmtree(triangle_dir, ignore_errors=True)
    print(f"Installing Triangle from {triangle_archive}...")
    shutil.unpack_archive(triangle_archive, CONFIG["directories"]["lib_dir"])
    run_command(f"make CC={CONFIG['compiler']['CC']}", cwd=triangle_dir)
    shutil.copy(os.path.join(triangle_dir, "triangle"), CONFIG["directories"]["bin_dir"])

def install_tetgen():
    """Install TetGen."""
    tetgen_archive = os.path.join(CONFIG["directories"]["lib_dir"], CONFIG["software"]["tetgen_archive"])
    tetgen_dir = os.path.join(CONFIG["directories"]["lib_dir"], CONFIG["software"]["tetgen_dir_name"])
    check_archive_exists(tetgen_archive)
    shutil.rmtree(tetgen_dir, ignore_errors=True)
    print(f"Installing TetGen from {tetgen_archive}...")
    shutil.unpack_archive(tetgen_archive, CONFIG["directories"]["lib_dir"])
    run_command(f"make CXX={CONFIG['compiler']['CXX']}", cwd=tetgen_dir)
    shutil.copy(os.path.join(tetgen_dir, "tetgen"), CONFIG["directories"]["bin_dir"])

def install_e4d_hr():
    """Compile and install E4D-HR."""
    petsc_dir = os.path.join(CONFIG["directories"]["lib_dir"], CONFIG["software"]["petsc_dir_name"])
    petsc_arch = "arch-linux-c-opt"
    e4d_path = os.path.join(CONFIG["directories"]["bin_dir"], "e4d")
    shutil.rmtree(e4d_path, ignore_errors=True)
    print("Compiling and installing E4D-HR...")
    env = os.environ.copy()
    env["PETSC_DIR"] = petsc_dir
    env["PETSC_ARCH"] = petsc_arch
    src_dir = os.path.join(CONFIG["directories"]["e4d_dir"], "src")
    run_command("make clean", cwd=src_dir, env=env)
    run_command(
        f"make FC='{CONFIG['compiler']['MPIFC']}' "
        f"FFLAGS='{CONFIG['compiler']['FOPTFLAGS']} -r8 -heap-arrays'",
        cwd=src_dir,
        env=env
    )
    shutil.copy(os.path.join(src_dir, "e4d"), e4d_path)

def menu():
    """Display menu and return selected option."""
    print("\nE4D-HR Installation Menu:")
    print("[1] Install All")
    print(" 2  Install PETSc")
    print(" 3  Install E4D-HR")
    print(" 4  Install Triangle")
    print(" 5  Install TetGen")
    choice = input("Select an option ([1]-5): ").strip()
    return choice if choice in ["2", "3", "4", "5"] else "1"

def main():
    print("E4D-HR Installation Script")
    if not os.getenv("MKLROOT"):
        print("Error: MKLROOT is not set. Please set the MKLROOT environment variable.")
        print("Try sourcing the Intel environment variables script, e.g., 'source <intel-compiler-dir>/setvars.sh'")
        return
    os.makedirs(CONFIG["directories"]["lib_dir"], exist_ok=True)
    os.makedirs(CONFIG["directories"]["bin_dir"], exist_ok=True)

    choice = menu()
    reinstall = prompt_reinstall()

    if choice == "1":
        if reinstall["PETSc"]:
            install_petsc()
        if reinstall["Triangle"]:
            install_triangle()
        if reinstall["TetGen"]:
            install_tetgen()
        if reinstall["E4D-HR"]:
            install_e4d_hr()
    elif choice == "2" and reinstall["PETSc"]:
        install_petsc()
    elif choice == "3" and reinstall["E4D-HR"]:
        install_e4d_hr()
    elif choice == "4" and reinstall["Triangle"]:
        install_triangle()
    elif choice == "5" and reinstall["TetGen"]:
        install_tetgen()
    print("Installation completed successfully.")

if __name__ == "__main__":
    main()
